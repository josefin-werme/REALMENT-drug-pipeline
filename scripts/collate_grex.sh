maf="maf005"
eqtls="z165"
tissue="brain_nucleus_accumbens"

prefix=collate_grex-$maf-$eqtls-$tissue

tmpdir="/gpfs/scratch1/nodespecific/tcn1309/josefin.tmp/collate/"
jobdir="$tmpdir/$tissue"
if [[ ! -d $tmpdir ]]; then mkdir $tmpdir; fi
if [[ ! -d $jobdir ]]; then mkdir $jobdir; fi

WT=8:00 # 10 min for subst nigra (8360174, 8360556), but 1h for cortex (8277150)

### 60k chunks, 30 in parallel; 30k arrays in total
ncores=16
method="R"

if [[ $tissue == "brain_substantia_nigra" ]]; then
	mem=28G
else
	mem=56G # cortex / whole_blood
fi

sbatch --mem $mem -t $WT:00 -J $prefix -o slurm.$prefix.%A --cpus-per-task=$ncores collate_grex.job $maf $eqtls $tissue $ncores $method $jobdir
