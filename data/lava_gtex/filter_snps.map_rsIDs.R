arg = commandArgs(T); tissue=arg[1]; lim = as.numeric(arg[2]) / 100

# tissue="whole_blood"; lim=1.65; c=22
library(data.table)

all_genes = list()

for (c in 1:22) {
	
	d = fread(paste0("all_snps/",tissue,"/",tissue,".eqtl.chr",c,".stats"), data.table=F)

	# get row numbers of all gene annotations (mark the start of the 1Mb window)
	genes = subset(d, gene!=".") # use the rownames of next gene to define end of current gene
	start_rows = as.numeric(rownames(genes)) # rows indicating the start of each gene
	n_genes = nrow(genes)

	genes_filt = list() 
	snp_to_gene = list()

	for (i in 1:n_genes) {
		start = as.numeric(start_rows[i]) # start index for eqtls for current gene
		if (i == n_genes) {
			stop = nrow(d)
		} else {
			stop = as.numeric(start_rows[i+1]) - 1 # stop index for eqtls for current gene
		}
		id = genes$gene[i] # gene id

		tmp = d[start:stop,] # subset whole data to current gene
		if (nrow(subset(tmp, gene!=".")) > 1) { stop("more than one gene in subset") } # safety check that theres only ever one gene listed in subsetted data

		tmp = subset(tmp, statistic < -lim | statistic > lim) # filter out variants that dont reach significance threshols
		# NOTE: this has to be done for each gene separately, since if I just subset the whole data set I might lose the SNP that indicates the gene ID

		if (nrow(tmp) > 0) {
			# store data in lava_gtex readable format
			genes_filt[[id]] = tmp
			genes_filt[[id]]$gene[1] = id # add gene ID back to first SNP (since the orig snp may have been filtered out)

			# also map every SNP ID to gene (used to determine which chunks to read in to LAVA)
			# I did this instead of mapping based on 1Mb from gene TSS as it seemed more error proof
			snp_to_gene[[id]] = tmp[,c(1,4)]
			snp_to_gene[[id]]$gene = id
		}
	}
	# store list of all genes containing at least one snp
	all_genes[[c]] = data.frame(chr = c, gene = names(genes_filt))

	# format lava_gtex files and snp-to-gene map
	filtered = do.call(rbind, genes_filt)
	snps_mapped = do.call(rbind, snp_to_gene)

	# save
	write.table(filtered, paste0("z",lim*100,"/",tissue,"/",tissue,".eqtl.chr",c,".stats"), row.names=F, quote=F)
	write.table(snps_mapped, paste0("z",lim*100,"/",tissue,"/",tissue,".eqtl.chr",c,".rsids"), row.names=F, quote=F)
}

# combine genes across chromosomes
all_genes = do.call(rbind, all_genes)

# order by chromosomal loccation
gencode = fread("gencode.v26.genes_ordered.gtf", data.table=F, header=F)
if (! all(all_genes$gene %in% gencode[,3])) { stop("Error: Not all genes in gencode file") }

# match to order in gencode file
ordered = all_genes[match(gencode[,3], all_genes$gene),]
ordered = ordered[!is.na(ordered[,1]),]

# do another check that all genes still exist
if (! all(ordered$gene %in% all_genes$gene)) { stop("Error: Some genes missing after ordering") }

# save
write.table(ordered, paste0("z",lim*100,"/",tissue,"/",tissue,".eqtl.genes"), row.names=F, quote=F)

