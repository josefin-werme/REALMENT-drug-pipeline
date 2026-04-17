# LAST EDITED: 31st Oct 2024

# subsets the full lincs data to named compounds (acc to hamonize_ids() function in local script ~/DrugSignatures/LINCS2020/filter_signatures.R), and exemplar signatures
# note that the lincs2.h5 is already subsetted to exemplar signatures (see local script ~/DrugSignatures/signatureSearchProcessData.R), but adding here to be xplicit

suppressPackageStartupMessages({
	library(ExperimentHub); library(SummarizedExperiment); library(HDF5Array); library(annotate); library(org.Hs.eg.db); library(data.table)
})

## Read in mappped exemplar signatures
exem.sigs.mapped = fread("exemplar_signatures_mappedIDs.txt", header=T, data.table=F)

## Read in lincs data
db_path = file.path("lincs2.h5")
se = SummarizedExperiment(HDF5Array(db_path, name="assay"))
rownames(se) = HDF5Array(db_path, name="rownames")
colnames(se) = HDF5Array(db_path, name="colnames")

## Create compound info
rows = strsplit(se@colData@rownames, "__")
coldat = data.frame(compound = sapply(rows, "[[", 1),
					cell = sapply(rows, "[[", 2),
					trt = sapply(rows, "[[", 3),
					row.names = se@colData@rownames)

## Re-process data with compound info
se = SummarizedExperiment(HDF5Array(db_path, name="assay"), colData = coldat)
rownames(se) = HDF5Array(db_path, name="rownames")
colnames(se) = HDF5Array(db_path, name="colnames")

## Extract signature dataset
sigdat = assay(se); dim(sigdat)

## Subset
sigdat.mapped = sigdat[,colnames(sigdat) %in% exem.sigs.mapped$compound_id]

## Save
writeHDF5Array(sigdat.mapped, filepath = "lincs2_mapped.h5", name = "lincs2_mapped", with.dimnames = TRUE)
#write.table(as.data.frame(sigdat.mapped), "lincs2_mapped.dat")
#fwrite(as.data.table(sigdat.mapped), "lincs2_mapped.dt")

## Test read in
# mat = HDF5Array(file = "lincs2_mapped.h5", name = "lincs2_mapped")

