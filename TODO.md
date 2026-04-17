# data
- add at least a subset of the data/lava_gtex so you can include the filter_snps.R script too
-- upload to zenodo and ADD LINK to data/lava_gtex/README.md

# misc script dir
- PRIO: remove redundant code in validation*.R scripts
- remove / replace readme
- remove the check_imputed / check_pred dirs ? (or keep for completion)
- where is predict_drugs.check_completed.R ?? 
- check the "collated="$grex/all_collated" # all_ind collated into single df (and class N subset)"
--- is this really everyone? if so is it per gene?
- where to put utils scripts? 
- function rerun() in impute_grex.job
- remove all the scripts/classN_*txt files? (these are just lists of phenotype IDs for the analyses; not integral and can change)
-- but if so I should probably make a not in scripts
- explain what defaultSubjectExclusions.txt and defaultVariantExclusions_rsids.txt are 
- include directory tree overview
- make sure output dirs are automatically created if not already in pred_drugs.sh and validate.sh

# signatures
- see todo in signatures dir
