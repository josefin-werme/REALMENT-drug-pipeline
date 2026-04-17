# ---- SPLITTING UKB DATA BY IND ---- #
1. Script directory: split_ukb

# ----- IMPUTATION ----- #

## Old strategy (group level)
1. impute_grex
2. collate_grex
3. predict_drugs - scores each drug based on deviation of target genes
4. validate - evaluates which drugs are significantly differnt in their scores based on certain case / control groups

## New strategy (v1; group level)
### old one relied on suboptimal drug signature data; now I'm instead using continuous gene expression signatures
0. data/lava_gtex/filter_snps.sh 	# filters SNPs for imputation based on Z score
1. impute_grex
2. compute_dev_grex_cases		# computes the difference between cases and controls in the form of Z scores (mean1 - mean2 / sd)
3. match_signatures			# matches the case-control differences in gene expression with signature data


