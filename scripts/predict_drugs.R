arg = commandArgs(T); infiles = arg[1]; outprefix = arg[2]; method = arg[3]; n.cores = as.numeric(arg[4])
# infiles='infiles.txt'; outprefix='maf005-z165-brain_cortex-zhang'; store_drugIDs = T; method="zhang"

suppressPackageStartupMessages({ library(bettermc); library(HDF5Array); library(data.table); library(annotate); library(org.Hs.eg.db); library(Hmisc) })

### Load drug signatures
#p0 = proc.time(); p0
sigdat = as.data.table(HDF5Array(file = "../lincs2_mapped.h5", name = "lincs2_mapped"), keep.rownames=T)

# subset signatures to only those of class N
classN.sigs = c("rn",read.table("../classN_signatures.txt",header=F)$V1) # read in class N signatures
sigdat = sigdat[,..classN.sigs] # subset signature data
num.sig.ids = 1:(ncol(sigdat)-1) # create numeric signature ids matching those of full signature data set

# get signature IDs and count
signatures = colnames(sigdat)[-1] # accounting for rownames
n.sigs = length(signatures)
print(paste0("n.sigs: ", n.sigs))

grex_files = read.table(infiles, header=F)$V1
n.ind = length(grex_files)

# subset sigdat
grex.genes = fread("grex.genes", header=F)
grex.genes = suppressMessages(mapIds(org.Hs.eg.db, keys=gsub("\\..*","",grex.genes), keytype="ENSEMBL", column="ENTREZID"))
#grex.genes = grex.genes[!is.na(grex.genes)]

shared.genes = intersect(sigdat$rn, grex.genes) # subsetting to shared genes
grex.keep = grex.genes[match(shared.genes, grex.genes)]

setkey(sigdat, rn)
sigdat = sigdat[shared.genes]
if (!all(sigdat$rn == grex.keep)) { stop("Error: subsetting of signature data gone wrong") }

progress = ceiling(quantile(1:n.ind, seq(.1,1,.1)))

## Get input files
null.out = mclapply(1:n.ind, mc.cores=n.cores, mc.retry=-1, function(i) { # mc.retry=-1 seems to not make a difference now as job just fails, but keeping anyway just in case
	if (i %in% progress) { print(paste0("* Progress: ",names(progress[which(progress==i)]))) }

	### Read in GREx for current individual
	infile = grex_files[i]
	grex = fread(infile, data.table=T)
	IID = grex$IID
	colnames(grex) = grex.genes

	### Convert to entrez IDs
	grex = grex[,..grex.keep] # remove genes which didnt convert
	if (!all(colnames(grex) == sigdat$rn)) { stop("Error: subsetting of grex gone wrong") }
	grex = unlist(grex)

	### Match grex with drug signatures
	if (method == "spearman") {	
		drugs = data.frame(sig = num.sig.ids, est = NA, p = NA)
	} else {
		drugs = data.frame(sig = num.sig.ids, C = NA, Cmax = NA, score = NA)
	}
	#progress = ceiling(quantile(1:n.sigs, seq(.1,1,.1)))

	#p0 = proc.time(); p0
	for (j in 1:n.sigs) {
		s = signatures[j]

		if (method == "spearman") {
			c = rcorr(x = sigdat[[s]], y = grex, type=method)
			drugs$est[j] = signif(c$r[1,2], 4)
			drugs$p[j] = signif(c$P[1,2], 4) 
		}
		if (method == "zhang") {
			# rank genes
			N = length(grex)
			order_G = order(abs(grex), decreasing=T)
			rank_G = (N - order_G + 1)
			s_rank_G = rank_G * sign(grex)

			order_S = order(abs(sigdat[[s]]), decreasing=T)
			rank_S = (N - order_S + 1)
			s_rank_S = rank_S * sign(sigdat[[s]])

			C = sum(s_rank_G * s_rank_S)
			Cmax = sum(rank_G * rank_G) # note, this is same as sum(rank_G[order(rank_G)] * rank_S[order(rank_S)])
			score = signif(C / Cmax, 6)

			drugs$C[j] = C; drugs$Cmax[j] = Cmax; drugs$score[j] = c
		}
	}
	#p0-proc.time()

#	if (infile=="chunk1.zgrex") { write.table(signatures, paste0("signatureIDs.",outprefix), col.names=F, row.names=F, quote=F) }

	# STORE OUTPUT
	# full results file without the signature IDs, since these are long and will be the same for all individuals 
	fwrite(drugs, paste0(IID,"-",outprefix,".drugs"), row.names=F, quote=F, sep=' ')

	# significant at p <.05 (incl the signature IDs)
#	fwrite(subset(drugs, p < .05), paste0(IID,"-",outprefix,".drugs_p05"), row.names=F, quote=F, sep=' ')
})


