source settings.sh
maf="maf005"
eqtls="z165"
tissue="brain_cortex" 
tissue="brain_substantia_nigra" # 
method="spearman"
pheno_item="d20003-ATC4N"

mem=28G
n_cores=16
WT=5:00

ROC=T	# set to "cases" or "controls" if run to determine which ones should be processed

enrich=F
classN_cases=T # only relevant for enrich; if F then all non_classN people will be analysed instead

tmpdir="/path/to/tmp/validate"

#for tissue in $(echo "brain_substantia_nigra" "brain_cortex" "brain_amygdala" "whole_blood"); do
for tissue in "brain_caudate"; do

if [[ $ROC == T ]]; then
#	atc_code_file="classN_ROC.txt"
	phenos=$(cat "classN_ROC.sub.txt")
else
	phenos=$("classN_phenos.txt")
fi

phenos=("N03AF")
#phenos=("notpure_code4")

for phen in ${phenos[@]}; do
#	if [[ $phen == "N06AA" ]]; then continue; fi

	printf "\nTissue: $tissue, ATC: $phen\n\n" #, cores: $n_cores
	prefix=val-$maf-$eqtls-$tissue-$method-$phen

	sbatch -t $WT:00 --mem $mem -J $prefix -o slurm.$prefix.%A --cpus-per-task=$n_cores \
	validate.job $maf $eqtls $tissue $method $n_cores $pheno_item $phen $tmpdir $enrich $ROC $classN_cases
done

done
