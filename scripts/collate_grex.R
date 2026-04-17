arg = commandArgs(T); prefix = arg[1]; suffix = arg[2]; outfile = arg[3]; genelist = arg[4]; n.cores = as.numeric(arg[5])
# prefix="zgrex_files/maf005-z165-whole_blood-"; suffix=".zgrex"; outfile="maf005-z165-whole_blood.zgrex_full"; genelist="maf005-z165-whole_blood.gene_list"; n.cores=15
library(data.table)
library(bettermc)

iids = fread("IIDs.txt",header=F,data.table=F)
genes = fread(genelist, data.table=F, header=F)$V1

grex_full = mclapply(1:length(genes), mc.cores = n.cores, mc.retry=-1, function(i) {
	infile = paste0(prefix, genes[i], suffix)
	return(fread(infile, data.table=F, header=F))
})
grex_full = cbind(iids, do.call(cbind, grex_full))
colnames(grex_full) = c("IID", genes)

fwrite(grex_full, outfile, row.names=F, col.names=T, sep=' ')
