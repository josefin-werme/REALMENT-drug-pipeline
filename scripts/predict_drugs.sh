maf="maf005" 
eqtls="z165" # SNP z-score filtering threshold for eqtl data
#tissue="brain_amygdala"
tissue="brain_caudate"
#method="zhang"
method="spearman"
echo "** Tissue: $tissue **"

WT=2:00
#mem=112G; n_cores=22; n_ind_per_job=320 # subst_nigra (class N cases, ALL signatures)
#mem=112G; n_cores=16; n_ind_per_job=320 # cortex (class N cases, ALL signatures)
mem=46G; n_cores=32; n_ind_per_job=320 # amygdala (classN signatures & cases)
mem=46G; n_cores=16; n_ind_per_job=320 # caudate (classN signatures & cases)
#mem=84G; n_cores=32; n_ind_per_job=320 # cortex (classN signatures/notN_cases)

classN_cases="TRUE" # whether only classN cases should be analysed or everyone else
#DEBUG="TRUE"; arr_debug=1

# ----------

source settings.sh
prefix=$maf-$eqtls-$tissue-$method
grxpref=$maf-$eqtls-$tissue

## Define common temp dir for prediction analyses
tmpdir="/path/to/your/temp"

# copy relevant files if they arent there already
cp $data/signatures/classN_signatures.txt $tmpdir
cp $scripts/predict_drugs.R $tmpdir
if [[ ! -f $tmpdir/lincs2_mapped.h5 ]]; then cp $data/signatures/lincs2_mapped.h5 $tmpdir; fi

if [[ $classN_cases == "TRUE" ]]; then
	if [[ ! -f $tmpdir/$grxpref.zgrex_classN ]]; then 
		cp $collated/$grxpref.zgrex_classN.gz $tmpdir # only class N individuals
		pigz -d $tmpdir/$grxpref.zgrex_classN.gz
	fi
	iid_file="$ukb_phenos/classN_meds_IIDs.txt"
else
	prefix=$prefix-not_classN
	iid_file="$ukb_phenos/not_classN_IIDs.txt"

	if [[ ! -f $tmpdir/$grxpref.zgrex_full ]]; then
		printf "\n >> classN_cases=FALSE; copying full GREx file to temp"
		cp $ukb_phenos/not_classN_IIDs.txt $tmpdir
		cp $collated/$grxpref.zgrex_full.gz $tmpdir # full collated grex dataset
		pigz -d $tmpdir/$grxpref.zgrex_full.gz
	fi
fi
n_ind_tot=$(nrow $iid_file)

## Make output dir
outdir=$drugs_pred/$prefix
if [[ ! -d $outdir ]]; then mkdir $outdir; fi 

# check if outfiles exist 
if [[ ! -d $scripts/check_pred ]]; then mkdir $scripts/check_pred; fi
cd $scripts/check_pred
outdir=$drugs_pred/$prefix
outfiles=($(ls $outdir/*.drugs.tar.gz 2>/dev/null))

# if no outfiles, submit all chunks
# otherwise, extract remaining chunks to submit
if [[ ${#outfiles[@]} == 0 ]]; then
	arr_run="1"
	printf "\n\n ** FIRST RUN: SUBMITTING ONLY ONE CHUNK ** \n\n"
else
	# get start/stop indexes from the outfiles
	ls $outdir/*.drugs.tar.gz | cut -f2 -d. | sort -n | sed 's/-/ /' > $prefix.fin
	
	# extract missing and save as a formatted string of array indices to submit
	Rscript predict_drugs.check_completed.R $prefix.fin $prefix $n_ind_per_job $iid_file    # outfile = "all_fin", "error", or string of remaining arrays

	# get the array indices from R output
	arr_run=$(cat $prefix.arr_run)
fi

if [[ $DEBUG == "TRUE" ]]; then 
	arr_run=$(echo $arr_run | cut -f1 -d'-')
	echo "** WARNING: ONLY ANALYSING $arr_run"
fi


## Check array indices, and submit
if [[ $arr_run == "error" ]]; then
	printf "\n\n ** ERROR: Problem extracting array indices ** \n\n"
	exit

elif [[ $arr_run == "all_fin" ]]; then
	echo "ALL finished !!!"
else
	cd $scripts
	sbatch -t $WT:00 -J predict-$prefix -o slurm.predict-$prefix.%A_%a --cpus-per-task=$n_cores --mem=$mem --array=$arr_run \
		predict_drugs.job $maf $eqtls $tissue $method $n_cores $n_ind_tot $n_ind_per_job $tmpdir $classN_cases
fi


