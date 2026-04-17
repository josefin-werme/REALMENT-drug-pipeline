
## LAST RAN ON 18th MAY 2024 (edits on 31st Oct are only to make comments)

# library(BiocParallel)
library(signatureSearch)
library(ExperimentHub)
#library(SummarizedExperiment); 
library(rhdf5)
library(HDF5Array)
library(data.table)


######     Lincs2     ######
setwd("~/surfdrive/POSTDOC/Genes_to_Drugs/DrugSignatures/LINCS2020/")

siginfo_beta <- fread("siginfo_beta.txt", data.table = F) # signature meta data
# exemplar <- siginfo_beta %>% filter(pert_type=="trt_cp" & is_exemplar_sig == 1) # subsetting to exemplar signatures treated with compounds

exemplar <- subset(siginfo_beta, pert_type=="trt_cp" & is_exemplar_sig == 1) # jw: had to use subset cause %>% wasnt working

### Keep only named compounds (i.e.remove all that only have a BRD ID since those are not useful for validation anyway)
# Code that produced this data can be found in filter_compounds() function in "validation_functions.R" script
comp_filtered = fread("compoundIDs_filtered.txt", data.table=F)
exemplar_filtered = exemplar[exemplar$pert_id %in% comp_filtered$pert_id,]

new_cid <- paste(exemplar_filtered$pert_id, exemplar_filtered$cell_iname, rep("trt_cp", length(exemplar_filtered$cmap_name)), sep="__") # creating IDs

gctx2h5("level5_beta_trt_cp_n720216x12328.gctx", cid=exemplar_filtered$sig_id, new_cid=new_cid,
        h5file="lincs2_filtered.h5", by_ncol=5000, overwrite=TRUE) # version 2 only to check warnings

# DBpath <- "lincs2.h5"
# sedb <- SummarizedExperiment(HDF5Array(DBpath, name="assay"))
# rownames(sedb) <- HDF5Array(DBpath, name="rownames")
# colnames(sedb) <- HDF5Array(DBpath, name="colnames")

# TODO: check if trt_cp is differential expression relative to controls
# Finally, we apply a moderated Z-scoring procedure to replicate samples of each experiment (generally 3 replicates are made) to compute a weighted average signature. For several types of genetic perturbagens (see pert_type table below; for example, where multiple hairpins were used to target the same gene) we further collapse level 5 data to create a single consensus gene signature (CGS) that represents that perturbation.

warnings()
print(warnings())
