# GTEx data
The full unfilterd gtex files are assumed to be in 'all_snps'. Note that these have been formatted specifically to work with the LAVA version used here
Due to file size limits these have been omitted; download link will be provided soon

# SNP filtering
To subset the GTEx files to only SNPs that pass a certain significance threshold, use filter_snps.sh and set the desired z-score value

# gencode file 
### used to extract chromosomal order of genes, so that genes in the from same/adjacent blocks will be analysed together
### this file is used to order the gene list created by filter_snps.job/R
awk '$3=="gene" {print $1,$4,$10}' gencode.v26.annotation.gtf | sed -s 's/"//g' | sed -s 's/;//' > gencode.v26.genes_ordered.gtf

##### orig file can be downloaded from https://www.gencodegenes.org/human/release_26.html


