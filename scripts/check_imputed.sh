## Checks whether there were any genes that did not finish during imputaiton, so that these can be re-imputed
 
function check_imputed () {
	# tissue="brain_cortex"; eqtls=z165
	debug=$1 # this works with current if statements regardless of what its set to

	pref=$tissue-$eqtls
	gene_file="$gtex/$eqtls/$tissue/*genes"
	outdir=$grex_ind/$tissue
	fin=$check_imputed/$pref.fin
	miss=$check_imputed/$pref.miss

	if [[ $debug ]]; then echo "Debug mode; $debug $pref $gene_file"; fi

	# Get finished
	{ ls $outdir/*zgrex.gz $outdir/*failed | sed "s/.*${tissue}-//" | cut -f1,2 -d'.' > $fin; } 2>/dev/null
	if [[ $debug ]]; then awk 'NR!=1' $gene_file > $fin; fi # triggers 'all'

	n_fin=$(nrow $fin)
	n_tot=$(($(nrow $gene_file) - 1)) 
	
	if [[ $debug ]]; then echo "n_fin $n_fin n_tot $n_tot"; fi

	if [[ $n_fin == 0 ]]; then 
		# If there are no finished genes, the gene file used for the imputation is just the full list of genes
		echo "none" 
	elif [[ $n_fin == $n_tot ]]; then
		# If all are finished
		echo "all"
	else
		# Get missing
		# (first column is excl since it is the 'gene' header)
		# awk 'NR==FNR {i[$1]; next} !($2 in i) {print $2}' $fin $gene_file | awk 'NR!=1' > $miss
		awk 'NR==FNR {i[$1]; next} !($2 in i)' $fin $gene_file > $miss
		if [[ $debug ]]; then head $miss > $miss; fi # triggers 'error'
		n_miss=$(($(nrow $miss) -1)) # -1 for the header

		# Check that total N rows of finished/missing add up to gene file,
		# and that non shared genes between missing and finished equal same as missing
		non_shared=$(($(comm -23 <(sort $miss) <(sort $fin) | nrow) -1))

		if [[ $((n_fin + n_miss)) != $n_tot  || $non_shared != $n_miss ]]; then 
			# If safety checks failed, return error
			echo "error"
		else
			# If safety checks passed and there is only a subset of genes to analyse:
			# .. the gene list for the analysis will be the $miss file
			echo "some"
		fi 
	fi
}
