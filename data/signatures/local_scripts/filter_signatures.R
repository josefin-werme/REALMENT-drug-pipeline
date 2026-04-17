
# --------------------------------------------------- #
# ----- SCRIPT TO FILTER DRUG DATA ADDITIONALLY ----- #
# --------------------------------------------------- #
# Last edited: 31st Oct 2024

# This script subsets the data to named compounds (acc to hamonize_ids() function), and exemplar signatures
# The initial processing/filtering script found in signatureSearchProcessData.R

# Create the filtering files locally, but RUN THE FINAL FILTERING ON THE CLUSTER (code at bottom) 

## Packages
suppressPackageStartupMessages({
	library(ExperimentHub); library(SummarizedExperiment); library(HDF5Array); library(annotate); library(org.Hs.eg.db); library(data.table)
})

## Directories
base.dir = "/path/to/your/base/dir"
lincs.dir = "/path/to/your/lincs/data"

# Read in ATC data and source drug ID harmonising function
source(file.path(base.dir, "validation_functions.R"))
read_atc()

#####
# -------------------------------------------------- #
# ------ FILTER COMPOUNDS TO ONLY NAMED DRUGS ------ #
# -------------------------------------------------- #
#####
# i.e. those that can be mapped using harmonize_ids() function
# this file can also be used to map the compounds to named drugs later

## Read in compound meta data
comp.all = fread(file.path(lincs.dir, "compoundinfo_beta.txt"), data.table = F)
comp.all = comp.all[,c("pert_id","cmap_name","compound_aliases","moa")]
dim(comp.all); head(comp.all)

## Harmonise drug IDs and subset to mapped drugs
univ.ids = harmonise_ids_vec(comp.all$cmap_name)
univ.alias = harmonise_ids_vec(comp.all$compound_aliases)

#-- where cmap conversion failed, use converted alias
sum(is.na(univ.ids))
na.univ.ids = is.na(univ.ids)
univ.ids[na.univ.ids] = univ.alias[na.univ.ids]
sum(is.na(univ.ids))

#-- set drug_ids to univ.ids
comp.all$drug_id = univ.ids
dim(comp.all); head(comp.all)

## Subset to only drugs with a mapped ID
comp.mapped = comp.all[!is.na(comp.all$drug_id),]
dim(comp.mapped); head(comp.mapped)

## Create final ID conversion file with only mapped drugs (and no duplicates)
comp.mapped = comp.mapped[,c('pert_id','drug_id','moa')]
comp.mapped = unique(comp.mapped)
dim(comp.mapped); head(comp.mapped)
write.table(comp.mapped, file.path(lincs.dir,"compoundIDs_mapped.txt"), row.names = F, quote = F, sep = '\t')


##### 
# ------------------------------------------------ #
# ---- Filter only on high quality signatures ---- #
# ------------------------------------------------ #
##### 
# NOTE: Filtering on the exemplar signatures only. The original downloaded and processed lincs2 data (read in below)
# was already filtered on these (see signatureSearchProcessData.R), 
# but still keeping this here for completion (and reassurance). 

siginfo.full = fread(file.path(lincs.dir, "siginfo_beta.txt"), data.table=F)

## Filtered data subsets
siginfo.exm = subset(siginfo.full, is_exemplar_sig == 1) # Exemplar signatures. 

## Set to chosen filter
siginfo = siginfo.exm

## Create compound ID equivalent to that in the lincs dataset
siginfo$compound_id = apply(siginfo[,c("pert_id","cell_iname","pert_type")], 1, paste, collapse="__")

## Create siginfo subset by filtering on named compounds (from above)
siginfo.mapped = subset(siginfo, pert_id %in% comp.mapped$pert_id)[,c("pert_id","cell_iname","pert_type","compound_id")]
dim(siginfo.mapped); head(siginfo.mapped)

## Add drug IDs to siginfo.mapped
siginfo.mapped$drug_id = comp.mapped[match(siginfo.mapped$pert_id, comp.mapped$pert_id),]$drug_id
head(siginfo.mapped)

## Save file (use to subset data on cluster)
write.table(siginfo.mapped, file.path(lincs.dir, "exemplar_signatures_mappedIDs.txt"), row.names=F, quote=F, sep = '\t')

##### 
# ------------------------- #
# ---- LOAD SIGNATURES ---- #
# ------------------------- #
##### 
# NOTE: these data sets were created using  signatureSearchProcessData.R

db_path = file.path(lincs.dir,"lincs2.h5") 

## Read in data
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

## Filter signatures
#--- check just the exemplar filter
sigdat.exm = sigdat[,colnames(sigdat) %in% siginfo$compound_id]
dim(sigdat.exm) # SAME as orig (i.e. reassurance that the lincs2.h5 data was already filtered on these)

#--- filter on both exemplar and compounds with mapped drug IDs 
sigdat.mapped = sigdat[,colnames(sigdat) %in% siginfo.mapped$compound_id]
dim(sigdat.mapped)

## Store (NOTE: do this only on cluster, code below)
# write.table(sigdat.hiq, file.path(lincs.dir,"lincsh5.comp-hiq.dat"))



##### 
# ------------------------------------------------ #
# ---- CODE TO FILTER THE DATA ON THE CLUSTER ---- #
# ------------------------------------------------ #
##### 
# NOTE: Doing the actual filtering on the cluster (rather than locally); this is the code
# Script saved on cluster as "filter_lincs_exemplar_mapped.R"

suppressPackageStartupMessages({
	library(ExperimentHub); library(SummarizedExperiment); library(HDF5Array); library(annotate); library(org.Hs.eg.db); library(data.table)
})

## Read in mappped exemplar signatures
sigs.mapped = fread("exemplar_signatures_mappedIDs.txt", header=T, data.table=F)

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
sigdat.mapped = sigdat[,colnames(sigdat) %in% sigs.mapped$compound_id]
sigdat.mapped

# writeHDF5Array(sigdat.mapped, filepath = "lincs2_mapped.h5", name = "lincs2_mapped", with.dimnames = T)
# write.table(sigdat.mapped, "lincs2_mapped.dat")

## NOTE: Script saved on cluster as "filter_lincs_exemplar_mapped.R"

# sm = sigdat.mapped[1:20,1:20]
# writeHDF5Array(sm, filepath = "lincs2_mapped.h5", name = "lincs2_mapped", with.dimnames = TRUE)
# mat = HDF5Array(file = "lincs2_mapped.h5", name = "lincs2_mapped")






























