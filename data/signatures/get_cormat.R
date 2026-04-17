# get covariance of signatures
library(data.table);library(HDF5Array); library(Hmisc)
sigdir = "/gpfs/work5/0/vusr0748/realment/data/signatures"
sigdat = as.data.table(HDF5Array(file = file.path(sigdir,"lincs2_mapped.h5"), name = "lincs2_mapped"), keep.rownames=F)
classN.sigs = read.table(file.path(sigdir,"classN_signatures.txt"), header=F)$V1 # read in class N signatures
sigdat = sigdat[,..classN.sigs] # subset signature data
c = rcorr(as.matrix(sigdat))
write.table(c$r, file.path(sigdir,"cormat_sigs.dat"), quote=F)

c = rcorr(t(as.matrix(sigdat)))
write.table(c$r, file.path(sigdir,"cormat_genes.dat"), quote=F)
