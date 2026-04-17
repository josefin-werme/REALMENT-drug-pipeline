eqtls="z165"
tissue="brain_putamen" 		# NOTE: last time there was 350 genes left to impute
tissue="brain_cerebellum" 		# NOTE: last time there was 350 genes left to impute
#tissue="pituitary" 		# NOTE: last time there was 350 genes left to impute
#tissue="brain_anterior_cingulate"
#tissue="brain_frontal_cortex"
#tissue="brain_hippocampus"
tissue="brain_hypothalamus"

#Wdebug="TRUE"; n_debug=10
#HIGH_MEM="TRUE"

# TODO: 'brain_anterior_cingulate/brain_anterior_cingulate_b24.eqtl.chr[CHR].stats'
#brain_cerebellum brain_frontal_cortex/brain_frontal_cortex_b9.eqtl.chr[CHR].stats
# TODO: why are not the correct n/o genes showing when submitting teh imputation???

npar=10
genes_per_job=50
mem="84G"
WT=6:00 
n_retry=0 


if [[ $debug == "TRUE" ]]; then 
	npar=$n_debug
	genes_per_job=$n_debug
fi

## ... high mem genes ... ##
if [[ $tissue == "whole_blood" ]]; then
	mem="224G" # max for the really big genes I saw was 64, but almost all were just under 60
	WT=10:00 # was 3 initially, doing
	npar=1 # 3
	genes_per_job=1 # 50
	n_retry=0 # 2
fi

# used for brain_cortex & brain_amygdala
if [[ $HIGH_MEM == "TRUE" ]]; then
        mem="100G" # max for the really big genes I saw was 64, but almost all were just under 60
        WT=8:00 # was 3 initially, doing
        npar=1 # 3
        genes_per_job=1 # 50
        n_retry=0 # 2
fi


# Check how many genes have already been imputed and adapt gene list accordingly
source settings.sh
imputed=$(check_imputed)

if [[ $imputed == "none" ]]; then
	gene_list="$gtex/$eqtls/$tissue/$tissue.eqtl.genes" # full gene list
elif [[ $imputed == "some" ]]; then
	gene_list="$check_imputed/$tissue-$eqtls.miss" # only missing
elif [[ $imputed == "all" ]]; then
	echo "All genes already imputed"; exit
elif [[ $imputed == "error" ]]; then
	"Error: Imputed genes cannot be determined"; exit
else
	"Error: Invalid return value from imputation check"
fi

# get total N genes 
n_genes=$(($(nrow $gene_list)-1)) # minus one for header
echo "N/o genes left to impute: $n_genes"

if [[ $debug == TRUE ]]; then
	n_genes=$n_debug
	echo "*** WARNING: DEBUG MODE, only $n_genes will be imputed ***"
fi

#for n_retry in 5; do
	#for npar in 3; do
		prefix=impute-${tissue}-${eqtls}-npar${npar}-retry${n_retry}
		sbatch --mem=$mem -t $WT:00 -J $prefix -o slurm.$prefix.%A_%a -n $((npar+1)) --array=1-$n_genes:$genes_per_job impute_grex.job $tissue $eqtls $npar $genes_per_job $n_retry $gene_list
	#done

#done
