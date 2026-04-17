# Script for converting SNP IDs

arg = commandArgs(T); bim_fname = arg[1]; ids_fname = arg[2]
# bim_fname = "chunk_206.bim"; ids_fname = "206.unique_to_rsid" #"ukb_unique_rsid_9.txt"

library(data.table)
ids = fread(ids_fname, header=F, data.table=F)
bim = fread(bim_fname, header=F, data.table=F)

### match IDs and check it is correct
matched.ids = ids[match(bim[,2],ids[,2]),]

print("***** Converting unique IDs to rsIDs - checking if done correctly *****")
if (!(all(matched.ids[,2] == bim[,2]))) { stop("Conversion of SNP IDs in bim file failed") }
# print(head(matched.ids[,1] == bim[,2]))

### replace orig ids with new ids in bim file
bim[,2] = matched.ids[,3]
#head(bim)
write.table(bim, paste0("rsids_",bim_fname), col.names=F, row.names=F, quote=F, sep='\t')

