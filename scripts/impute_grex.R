arg = commandArgs(T); tissue = arg[1]; genes = unlist(strsplit(arg[2], ";")); chrs = unlist(strsplit(arg[3], ";")); outname = arg[4]; ncores = as.numeric(arg[5]); retry.flag = as.numeric(arg[6])
# tissue='brain_cortex'; ncores=5; genes="ENSG00000263513.5"; chrs=1; outname="test.debug"

# genes="ENSG00000263513.5;ENSG00000162782.15;ENSG00000143340.6;ENSG00000135837.15;ENSG00000206527.9;ENSG00000065534.18;ENSG00000175455.14;ENSG00000065371.17;ENSG00000160145.15;ENSG00000151552.11;"; genes = unlist(strsplit(genes, ";")); chrs=1; tissue='brain_cortex'; outname="debugtest"; ncores=1; retry.flag=0

library(bettermc)
library(data.table)
# remotes::install_local("/home/josefin/programs/LAVA_batch_locus/", force=T)
pkgload::load_all("/home/josefin/programs/LAVA_batch_locus/")


## --- INPUT PROCESSING --- ##
#
## Load and validate data files files for specified chromosome (creates a DataSet object)
# ... input.dir is appended as prefix to file names specified in info file
data = process.input(input.info="gtex_v8.info", sample.overlap.file=NULL, 
			ref.prefix="chunk_{BLOCK}_processed", chromosomes=chrs, 
			input.dir=tissue, phenotypes=tissue)

## Create DataInterface object for all phenotypes in DataSet object, and
# ... create internal annotation object based on the 'gene' unit from the input files
input = data$get.interface(tissue, minimum.snps = 1, minimum.components = 1)$create.annotation("gene")
# NOTE: locus.miminum.lower.bound has been changed in R/globals.R to 1 instead of 5 since I would like to be able to impute data for genes with only a single SNP


#### GENE LIST ####
# ... since we only want to analyse the genes passed from the job file, we can't just iterate through the genes in the DataInterface object.
# ... there might also be some genes which for some reason didnt pass QC, so these should be moved from our gene list as well
processed.genes = input$get.annotation()$get.names()
processed.genes = processed.genes[processed.genes %in% genes]
n.genes = length(processed.genes)

# ... save empty output file for genes which werent processed
if (length(genes) != n.genes) {
	for (curr.gene in genes[! genes %in% processed.genes]) {
		write.table(NA, paste0("out.tmp/",outname,"-",curr.gene,".procfailed"), row.names=F, col.names=F)
	}
}

# ... if no genes remain, exit
if (n.genes == 0) {
	stop("None of the current genes passed QC")
}


# output data frames
grex = list(); locus.info = list()

# progress printing
progress = round(quantile(1:n.genes, seq(.05,1,.05)))
print_parallel <- function(...){
	system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
}

n.dec = 4 # n/o decimals to store output

# loop over genes
# preschedule=F ensures the virtual memory doesnt keep scaling with each itr
null.out = mclapply(1:n.genes, mc.retry = retry.flag, mc.silent = F, mc.messages = "m_ignore", mc.cores = ncores, mc.preschedule=F, function(i) {
#  for (i in 2:n.genes) { 
#print(""); print(paste0("**** CURRENT ITERATION = ",i," ****")); print("")
	if (i %in% progress) { print_parallel(paste0("Progress: ",names(progress)[progress %in% i])) }
	
	curr.gene = processed.genes[i]
 
	input$set.locus(curr.gene)
	processed.locus = catchErrors(input$process(compute.marginal=F)$process(), all.error=T)
	gc() 

	if (!is.null(processed.locus)) {
	    if (processed.locus$no.pheno("available") > 0) { 

		phen.id = processed.locus$phenotypes()

		# store locus info
		locus.info = data.frame(gene = curr.gene)
		locus.info$omega = signif(processed.locus$get.estimates()[[phen.id]]$get("omega"), n.dec)
		locus.info$sigma = signif(processed.locus$get.estimates()[[phen.id]]$get("sigma"), n.dec)
		locus.info$n.pcs = processed.locus$no.components()
		locus.info$no.snps = processed.locus$no.snps()
		fwrite(locus.info, paste0("out.tmp/",outname,"-",curr.gene,".locus"), sep="\t")

		# compute grex
		w = processed.locus$get.ld()$get.components()
		delta = processed.locus$get.estimates()[[phen.id]]$get("delta")
		grex = signif(w %*% delta, n.dec)
		rm (w, delta, phen.id, processed.locus); gc()

		# compute z grex
		z.grex = signif((grex - mean(grex)) / sd(grex), n.dec)

		fwrite(grex, paste0("out.tmp/",outname,"-",curr.gene,".grex"), col.names=F, row.names=F)
		fwrite(z.grex, paste0("out.tmp/",outname,"-",curr.gene,".zgrex"), col.names=F, row.names=F)

		rm(grex, z.grex); gc()

	    } else {
		rm(processed.locus); gc() # dont know if necessary but just in case
	        write.table(NA, paste0("out.tmp/",outname,"-",curr.gene,".failed"), row.names=F, col.names=F)
		}
	} else {
	    rm(processed.locus); gc()
	    write.table(NA, paste0("out.tmp/",outname,"-",curr.gene,".failed"), row.names=F, col.names=F)
	}
	return(NULL)
}

)




