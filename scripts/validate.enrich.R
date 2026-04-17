arg = commandArgs(T); input_dir = arg[1]; inprefix = arg[2]; outdir = arg[3]; code = arg[4]; case_file = arg[5]; n.cores = as.numeric(arg[6]); only.first.test = arg[7]
# code="N06BA"; tissue="brain_cortex"; input_dir="drugs"; inprefix=paste0("maf005-z165-",tissue,"-spearman"); case_file = paste0("../classN_cases/",code,".cases"); n.cores = 15; outdir = paste0("enrich_atc/",code)

# for cortex, check "N06AA" and "N05AB"
library(data.table); library(bettermc)

### Signature IDs and annotation
sigs = read.table("../signatureIDs_mapped.txt", header=F); colnames(sigs) = "id"
sig.info = fread("../exemplar_signatures_mappedIDs.txt", data.table=F)
classN.sigs = fread("../classN_signatures.txt",data.table=F,header=F)$V1 # NOTE: subsetting to classN sigs since only these were analysed since whole_blood (cortex and nigra have all)
sigs = data.frame(id = sigs[match(classN.sigs, sigs$id),])
sigs$num = 1:nrow(sigs) # numeric id equivalent to that in results
sigs$drug = sig.info[match(sigs$id, sig.info$compound_id),]$drug # drug id

### ATC classes
atc = fread("../atc_processed.txt", data.table=F)
atc = subset(atc, drug %in% sigs$drug) # subset to drugs in signature data to avoid analysing atc classes without any drugs 
atc = subset(atc, code1 == "N") # subset of class N # NOTE: subsetting to classN

level = nchar(code)-1
code.level = paste0("code",level) # NOTE: this wont work with class1, but maybe wont be relevant?
class = unique(atc[atc[[code.level]] == code, code.level])
if (length(class)!=1) { stop("ATC file incorrectly processed") }

atc.codes = unique(atc[[code.level]]) # all unique atc codes to evaluate enrichment
n.atc = length(atc.codes)

# extract list of drugs for each atc class to be used in the enrichment 
atc.drugs = mclapply(atc.codes, mc.cores=n.cores, function(j) {
	unique(atc[atc[[code.level]] == j,]$drug)
})
names(atc.drugs) = atc.codes

### Case IIDs and relevant infiles
cases = read.table(case_file, header=F)$V1 # IIDs for cases
infiles = sapply(cases, function(x) paste0(x,"-",inprefix,".drugs")) # make vector of relevant infiles from case IIDs

## Moa
moa = fread("../compoundIDs_mapped.txt", data.table=T)[,-1]

## Signatures/drugs listed in atc file
sigs.enrich = subset(sigs, drug %in% atc$drug)$num 
drugs.enrich = sigs$drug[match(sigs.enrich, sigs$num)] # drug ids for signatures that overlap with drugs in atc file

### Read in drug prediction results
print("* Reading in data ...")
if (code=="N06A") { n.cores = 14 }

pred = new.env()
#drugs.pred = mclapply(1:length(infiles), mc.cores = n.cores, function (i) {
pred$sigs = mclapply(1:length(infiles), mc.cores = n.cores, function (i) {
	dat = try(fread(file.path(input_dir,infiles[i]), data.table=T), silent=T)
	if (class(dat)[1] != "try-error") { 
        	dat = subset(dat, sig %in% sigs.enrich) # subset predicted signatures to those overlapping with the ATC file
        	dat$drug = drugs.enrich # add drug id column
		setkey(dat, p)
		return(dat) 
	}
})
# remove null entires (these are from people who were in total ukb data set but later filtered)
null.dat = sapply(1:length(pred$sigs), function(x) is.null(pred$sigs[[x]]))
pred$sigs = pred$sigs[!null.dat]
n.pred = length(pred$sigs) # same as N but keeping now anyway cause I have used in in old code
cases = cases[!null.dat]
N = length(cases)
gc()

# MEMORY: 7.6 Gb
# ----------------------------------------------------------------------------------------------------------- #
#     Check how often at least one drug from relevant category is prioritised (only looking at category N)    #
# ----------------------------------------------------------------------------------------------------------- #
# check what is the top drug of any N
# rank of relevant drug
get.ranks=F
if (get.ranks==T) {
	ranks = mclapply(atc.codes, mc.cores = n.cores, mc.retry = -1, function(code.j) { 
		# for every category, get rank?
		x = list()
		x$top = rep(NA, length(N))
		x$median = rep(NA, length(N))

		for (i in 1:N) {
			idx = which(pred$sigs[[i]][drug %in% unlist(atc.drugs)]$drug %in% atc.drugs[[code.j]])
			#x$idx[[i]] = idx
			x$top[i] = min(idx)
			x$median[i] = median(idx)
		}
		return(x)
	})
	names(ranks) = atc.codes

	summary = list()
	vars = c('top','median')
	for (v in vars) {
		summary[[v]] = apply(sapply(ranks, "[[", v), 2, mean)
		summary[[v]] = signif(summary[[v]][order(summary[[v]])],2)
		write.table(summary[[v]], paste0("ranks/",code,".",v), row.names=T, col.names=F, quote=F)
	}
}


# ------------------------------------------------ #
#    Evaluate enrichment of significant results    #
# ------------------------------------------------ #
print("* Evaluating enrichment for every individual separately ...")
p.thresh.str = "5e-2"
p.thresh = as.numeric(p.thresh.str)

## Enrichment for every individual separately
e.atc = list()
progress = round(quantile(1:n.atc, seq(.1,1,.1)))

if (code=="N06A") { n.cores = 5 }

#if (code !="N03A" & (N < 2500 | code.level != "code4")) { print ("*** WARNING *** only analysing individual level enrichment for N < 2500") # NOTE: <<<<<< DEBUG <<<<<
run.ind = T
if (run.ind) {
	p0 = proc.time()
	for(j in 1:n.atc) {
		code.j = atc.codes[j]
		drugs.j = atc.drugs[[code.j]]

		# uncomment this to read in existing e.atc data
		#	for (i in 1:n.pred) { e.atc[[j]] = fread(paste0(outdir,"/",code,"_",code.j,"-",inprefix,".enrich_",p.thresh.str)) }}

		e.atc[[code.j]] = do.call(rbind, mclapply(1:n.pred, mc.cores=n.cores, mc.retry=-1, function(i) {
			# check if there are any significant drugs for this atc code
			# evaluate enrichment only for individuals with at least one sig drug
			if (any(pred$sigs[[i]][p < p.thresh]$drug %in% drugs.j)) {
				sig = subset(pred$sigs[[i]], p < p.thresh)$drug %in% drugs.j
				nsig = subset(pred$sigs[[i]], p > p.thresh)$drug %in% drugs.j
				fishdat = data.frame(
					"sig" = c(sum(sig == T), sum(sig == F)),
					"not_sig" = c(sum(nsig == T), sum(nsig == F)),
					row.names = c("related", "not_related"),
					stringsAsFactors = FALSE
				)
				fish = try(fisher.test(fishdat, alternative = "greater"), silent=T)
				data.table(
					iid = cases[i],
					or = round(as.numeric(fish$estimate), 2),
					p = signif(fish$p.value, 2))
			} else {
				# otherwise just set to 0,1
				data.table(
					iid = cases[i],
					or = 0,
					p = 1)
			}
		}))
		fwrite(e.atc[[j]], paste0(outdir,"/",code,"_",code.j,"-",inprefix,".enrich_",p.thresh.str), sep="\t")
		if (j %in% progress) { print(paste0("* Progress: ",names(which(progress==j)))) }
	}
	proc.time() - p0
	invisible(gc())


	# get top atc classes based on pvals
	for (sig.thresh in c("5e-2","5e-3","5e-4")) {

		n.sig = data.frame(atc = names(e.atc), n.sig = sapply(1:length(e.atc), function(i) sum(e.atc[[i]]$p < as.numeric(sig.thresh))))
		#setkey(n.sig, n.sig) # order by p
		n.sig = n.sig[base::order(n.sig$n.sig, decreasing=T),]

		n.sig = cbind(n.sig, do.call(rbind,
			lapply(n.sig$atc, function(x) {
				atc.sub = atc[atc[[code.level]] == x,][,c(paste0("class",1:level),"drug")] # only include classes up until current code level
				classes = unique(atc.sub[,-ncol(atc.sub)])
				if (nrow(classes) > 1) { print("ERROR") }
				drugs = unique(atc.sub$drug)
				classes$drugs = paste(drugs, collapse="; ")
				moa.x = unique(moa$moa[match(drugs, moa$drug_id)])
				moa.x = moa.x[moa.x!=""]
				classes$moa = paste(moa.x, collapse="; ")
				return(classes)
			})
		))
		n.sig = subset(n.sig, n.sig != 0)
		if (nrow(n.sig) > 0) { write.table(n.sig, paste0("enrich_atc/",code,"-",inprefix,".","nsig_",sig.thresh), sep="\t", row.names=F, quote=F) }
	}
} 


## Enrichment of average results across individuals
# OPTION 1: extract sigs above/below p.thresh from all ind
# subset to those in atc, and map to drugs
# iterate over atc and evaluate enrichment
print("* Evaluating enrichment for all individuals jointly ...")

# extract sigs above/below p.thresh from all ind 
pred$sig = unlist(sapply(1:n.pred, function(i) pred$sigs[[i]][p < p.thresh]$sig))
pred$nsig = unlist(sapply(1:n.pred, function(i) pred$sigs[[i]][p > p.thresh]$sig))
# subset to those in atc
pred$sig = pred$sig[pred$sig %in% sigs.enrich]
pred$nsig = pred$nsig[pred$nsig %in% sigs.enrich]
# map to drugs
pred$sig = sigs$drug[match(pred$sig, sigs$num)]
pred$nsig = sigs$drug[match(pred$nsig, sigs$num)]

# MEMORY: 14.3 Gb
#sapply(ls(), function(x) format(object.size(get(x)), units="auto"))
#sapply(ls(pred), function(x) format(object.size(get(x, envir = pred)), units="auto"))





#### Enrichment

# --- OPTIMISATION ---
# Below enrichment runs out of memory when N >= 2026 - 3762 somewhere
# TODO: optimise by using numeric drug ids & set cores based on N

## Rationale:
# 2026 * 16 = 32416, therefore
# 32416 / N = n.cores

n.cores = round(32416 / N)
n.cores = n.cores - 2
n.cores = min(16, n.cores) # set to a max of 16
n.cores = max(1, n.cores) # and min of 1

e = do.call(rbind, mclapply(1:n.atc, mc.cores=n.cores, mc.retry=-1, function(j) {
	code.j = atc.codes[j]
        drugs.j = atc.drugs[[code.j]] # extract drugs for current atc class

        sig = pred$sig %in% drugs.j
        nsig = pred$nsig %in% drugs.j
	fishdat = data.frame(
                "sig" = c(sum(sig == T), sum(sig == F)),
                "not_sig" = c(sum(nsig == T), sum(nsig == F)),
                row.names = c("related", "not_related"),
                stringsAsFactors = FALSE
        )
        fish = try(fisher.test(fishdat, alternative = "greater"), silent=T)

	if (j %in% progress) { print(paste0("* Progress: ",names(progress[which(progress==j)]))) }
	data.frame(
		atc = code.j,
		or = round(as.numeric(fish$estimate), 2),
                p = signif(fish$p.value, 2))
}))
gc()

# Format enrichment results
e = e[order(e$p),]
e = cbind(e, do.call(rbind, 
	lapply(e$atc, function(x) {
		atc.sub = atc[atc[[code.level]] == x,][,c(paste0("class",1:level),"drug")] # only include classes up until current code level
		classes = unique(atc.sub[,-ncol(atc.sub)])
		if (nrow(classes) > 1) { print("ERROR") }
		drugs = unique(atc.sub$drug)
		classes$drugs = paste(drugs, collapse="; ")
		moa.x = unique(moa$moa[match(drugs, moa$drug_id)])
		moa.x = moa.x[moa.x!=""]
		classes$moa = paste(moa.x, collapse="; ")
		return(classes)
	})
))

write.table(e, paste0(code,"-",inprefix,".enrich-all"), row.names=F, quote=F, sep='\t')

#rownames(e) = 1:nrow(e)
#subset(e, atc == code)
#head(subset(e, class1 == "Nervous System"))

# TODO: ?? enrichment of moa based on e ??


# OPTION 2: rank drugs based on pvalues and correlation (only for ties)
#pvals = lapply(pred$sigs, "[[", "p")
#pvals = sapply(1:length(pvals[[1]]), function(row) mean(sapply(pvals, "[", row)))

#pred$avg = data.frame(sig = sigs$id, drug = sigs$drug, p = pvals)
#pred$avg = pred$avg[order(pred$avg$p),]
#write.table(pred$avg, paste0(outprefix,".sigs-avg-p"), row.names=F, quote=F)

enrich.env = function(env, df, idx, p.thresh, cores) {
        out = do.call(rbind,
                mclapply(atc.codes, mc.cores=cores, retry=-1, function(j) { # j = "N06AA"

                        drugs.j = atc.drugs[[code.j]] # extract drugs for current atc class
                        sig = subset(env[[df]][[idx]], p < p.thresh)$drug %in% drugs.j
                        nsig = subset(env[[df]][[idx]], p > p.thresh)$drug %in% drugs.j

                        fishdat = data.frame(
                                "sig" = c(sum(sig == T), sum(sig == F)),
                                "not_sig" = c(sum(nsig == T), sum(nsig == F)),
                                row.names = c("related", "not_related"),
                                stringsAsFactors = FALSE
                        )
                        fish = try(fisher.test(fishdat, alternative = "greater"), silent=T)
                        data.frame(
                                atc = code.j,
                                or = round(as.numeric(fish$estimate), 2),
                                p = signif(fish$p.value, 2))
                })
        )
	return(out)
}

enrich.df = function(df, p.thresh, cores) {
        out = do.call(rbind, 
		mclapply(atc.codes, mc.cores=cores, retry=-1, function(j) { # j = "N06AA"

			drugs.j = atc.drugs[[code.j]] # extract drugs for current atc class
			sig = subset(df, p < p.thresh)$drug %in% drugs.j
			nsig = subset(df, p > p.thresh)$drug %in% drugs.j

			fishdat = data.frame(
				"sig" = c(sum(sig == T), sum(sig == F)),
				"not_sig" = c(sum(nsig == T), sum(nsig == F)),
				row.names = c("related", "not_related"),
				stringsAsFactors = FALSE
			)
			fish = try(fisher.test(fishdat, alternative = "greater"), silent=T)
			data.frame(
				atc = code.j,
				or = round(as.numeric(fish$estimate), 2),
				p = signif(fish$p.value, 2))
		})
	)
	return(out)
}







