# Directories
ukb_chunks="/path/to/UKBdata"
project_dir="/path/to/project_dir"
scripts="$project_dir/scripts"
data="$project_dir/data" # all non genotype data
gtex="$data/lava_gtex" # eqtls processed for new LAVA version
grex="$project_dir/grex_imputed"
grex_ind="$grex/all_ind" # per gene GREx computed for all individuals in data set
collated="$grex/all_collated" # all_ind collated into single df (and class N subset)
drugs_pred="$project_dir/drugs_pred"
drugs=$drugs_pred
ukb_phenos="$data/ukb_phenotypes"
sigs="$data/signatures"
validation="$project_dir/validation"

# Utils
utils="$HOME/utils/scripts"
source $utils/check_disk_quota.sh
source $utils/count_rows.sh
source $utils/count_files.sh
check_imputed="$scripts/check_imputed"
source $scripts/check_imputed.sh

# Total sample size variable
ukb_n=387419

