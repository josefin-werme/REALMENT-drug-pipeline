
### Directory tree structure
The scripts assume the following directory structure
Note that all key paths are specified in scripts/settings.sh file for convenice. As long as you set the root dir the rest will be inferred. 

.
├── data
│   ├── lava_gtex
│   │   ├── all_snps
│   │   ├── z165
│   │   └── z196
│   ├── signatures
│   └── ukb_phenotypes
├── drugs_pred
│   ├── maf005-z165-brain_amygdala-spearman  # output file directories named acc to: ${maf}-${eqtl_snp_z_thresh}-${tissue}
│   ├── maf005-z165-brain_cortex-spearman
│   └── ...etc
├── grex_imputed
│   ├── all_collated
│   ├── all_ind
├── scripts
│   ├── check_imputed
│   └── check_pred
└── validation
    ├── maf005-z165-brain_amygdala-spearman
    │   ├── enrich_atc
    │   └── ROC
    ├── maf005-z165-brain_cortex-spearman
    │   ├── enrich_atc
    │   └── ROC
    └── ...etc
