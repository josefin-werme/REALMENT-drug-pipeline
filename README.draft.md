# Project summary
This repository contains the scripts and key resources used for the development of a drug prediction pipeline as part of REALMENT [..grant details..].

This pipeline relied on genome-wide genetic data & eQTL data from GTEx to impute genetically regulated gene expression across different tissues, matched with drug-expression signatures from [...] in order to predict individual level drug response.

A key part of this project was to determine the impact of [parameter choices/setup/datasets] to optimise prediction accuracy.

This included
- Selection of SNPs
- Selection of tissues
- Method used for signature matching
- [Drug signature dataset?]

# Repository overview
### programs
#### |- LAVA_batch_locus
Custom implementation of one of the earlier versions of the LAVA R package [.. details: TWAS, batch locus processing].
This package was used to impute the genetically regulated gene expression, based on tissue specific eQTLs from GTEx.
#### TODO:  >> just put programs in project dir instead of separate? << 

### project_dir
#### |- data
[key files]
#### |--- lava_gtex
Processed GTEx data in the appropriate format for the custom LAVA implementation
[different levels of SNP filtering]
#### |--- signatures

# Workflow
1. Filter eQTLs based on significance
  - [relevant scripts]
3. Impute gene expression for all individuals in the data set
  - [relevant scripts]
4. Predict individual level drug response via signature matching, separately for each available tissue. This was done using X differnet prediction methods: [...].
  - [relevant scripts]
5. Evaluate predictive performance
  - [relevant scripts]



