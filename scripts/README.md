# Pipeline Overview

Scripts implement a 4-step pipeline to impute genetically regulated gene expression (GREx) from UK Biobank genotypes, then use those imputed values to predict drug responses and validate predictions against true medication prescription info.

**Shared settings:** All scripts source `settings.sh`, which defines paths (`$grex`, `$collated`, `$drugs_pred`, `$sigs`, `$ukb_phenos`, etc.) and utility functions (`check_disc_quota`, `nrow`, `nfiles`).

---

## Step 1: Impute GREx — `impute_grex.sh` → `impute_grex.job` → `impute_grex.R`

**Goal:** For a given tissue, impute the genetically regulated component of gene expression (GREx) for every individual in the UKB cohort, for every gene with available GTEx eQTL data.

### `impute_grex.sh`
- Sets the target tissue, eQTL filtering threshold (`$eqtls`), and SLURM resource params (memory, walltime, parallelism).
- Checks which genes have already been imputed (via `check_imputed`): uses the full gene list, or only missing genes if a run was partially completed.
- Submits `impute_grex.job` as a SLURM array — one array task per gene (or batch of `genes_per_job` genes).

### `impute_grex.job`
- Receives: tissue, eqtls threshold, parallelism, gene list, retry flag.
- Creates a scratch temp dir; cleans up on exit.
- Checks disk quota before proceeding.
- Copies to temp: GTEx eQTL files for the tissue, UKB SNP↔rsID↔chunk mapping, exclusion lists, and R scripts.
- For each SLURM array task (a batch of genes):
  1. Identifies the relevant genotype chunks for the current genes' SNPs.
  2. For each chunk: copies plink files, converts UKB-unique IDs to rsIDs (`unique_to_rsid.R`), runs `plink2` to extract relevant SNPs (MAF ≥ 0.5%, applies subject/variant exclusions).
  3. Runs `impute_grex.R` in parallel across genes.
  4. If any genes have missing output, reruns with `mc_retry=-1` (progressively reduces cores on failure).
- Copies output files (`.grex.gz`, `.zgrex.gz`, `.locus`, `.failed`) to `$grex/all_ind/$tissue/`.
- On the first array task: saves the IID list for later collation.

### `impute_grex.R`
- Uses the custom `LAVA_batch_locus` package to impute GREx.
- Loads and validates genotype + eQTL data for the relevant chromosomes.
- For each gene (parallel via `bettermc::mclapply`):
  1. Sets the locus to the current gene, processes it (LD decomposition, effect estimation).
  2. Saves locus-level statistics (omega, sigma, n PCs, n SNPs) to a `.locus` file.
  3. Computes GREx: `w %*% delta` (LD components × effect estimates).
  4. Computes z-GREx: `(GREx − mean) / sd`.
  5. Writes `.grex` and `.zgrex` files; writes `.failed` if processing fails.

**Output per gene per individual:** raw GREx and z-scored GREx values.

---

## Step 2: Collate GREx — `collate_grex.sh` → `collate_grex.job` → `collate_grex.R`

**Goal:** Combine the per-gene, per-individual GREx files from Step 1 into a single wide-format matrix, then subset to class N medication users.

### `collate_grex.sh`
- Sets tissue, maf, eqtls; creates temp dir.
- Submits a single (non-array) SLURM job via `collate_grex.job`.

### `collate_grex.job`
- Copies all `.zgrex.gz` files for the tissue from `$grex/all_ind/$tissue/` to temp, decompresses.
- Creates a gene list from the filenames.
- Calls `collate_grex.R` to combine all per-gene zgrex files into one matrix (individuals × genes), with IIDs as the first column and gene names as column headers.
- Subsets the full matrix to class N medication cases (individuals in `classN_meds_IIDs.txt`) using `awk`.
- Compresses both the full and classN matrices with `pigz` and copies to `$collated/`.

### `collate_grex.R`
- Reads all per-gene `.zgrex` files in parallel (`bettermc::mclapply`).
- Column-binds them with IIDs → produces a single `individuals × genes` data frame.
- Writes to a space-delimited file.

**Output:** Two files in `$collated/`:
- `$prefix.zgrex_full.gz` — all individuals
- `$prefix.zgrex_classN.gz` — class N medication users only

---

## Step 3: Predict Drug Responses — `predict_drugs.sh` → `predict_drugs.job` → `predict_drugs.R`

**Goal:** For each individual, evaluate their GREx profile against LINCS2 drug-induced transcriptomic signatures to generate a ranked list of drugs.

### `predict_drugs.sh`
- Sets tissue, method (`spearman` or `zhang`), resource params, and `classN_cases` flag.
- Copies the relevant collated GREx file (classN or full) to a shared tmpdir and decompresses it.
- Copies LINCS2 signature file (`lincs2_mapped.h5`) and `predict_drugs.R` to tmpdir.
- Checks for existing output files and identifies which individual-chunks still need processing (`predict_drugs.check_completed.R`).
- Submits `predict_drugs.job` as a SLURM array, with each task processing `n_ind_per_job` (default 320) individuals.

### `predict_drugs.job`
- Receives: maf, eqtls, tissue, method, n_cores, n_ind_tot, n_ind_per_job, tmpdir, classN_cases.
- Computes the start/stop individual index for this array task.
- For each individual in the chunk: extracts their single row from the shared GREx matrix → saves as `chunk$J.zgrex`.
  - ClassN mode: subsets by row index.
  - Not-classN mode: subsets by IID match.
- Runs `predict_drugs.R` via `srun`.
- Verifies all output files exist, tars them (pigz-compressed) and copies to `$drugs_pred/$outname/`.

### `predict_drugs.R`
- Loads LINCS2 drug signatures from HDF5, subsets to class N signatures.
- Maps GREx gene IDs (Ensembl → Entrez), finds overlap with signature data.
- For each individual (parallel):
  1. Reads their `chunk$J.zgrex` file.
  2. Subsets to shared genes.
  3. For each drug signature, computes similarity:
     - **Spearman**: Spearman correlation (r, p-value) between individual GREx and drug signature.
     - **Zhang**: rank-based connectivity score (C/Cmax).
  4. Writes a per-individual `.drugs` file: one row per drug signature, with correlation estimate and p-value.

**Output:** Per-individual `.drugs` files, tarred by chunk, stored in `$drugs_pred/$prefix/`.

---

## Step 4: Validate — `validate.sh` → `validate.job` → `validate.roc.R` / `validate.enrich.R`

**Goal:** Evaluate whether drug predictions are enriched for the drugs individuals actually take (known medication data from UKB).

### `validate.sh`
- Sets tissue, method, ATC phenotype codes, resource params, and flags: `ROC=T`, `enrich=F`.
- Loops over tissues and ATC codes; submits one `validate.job` per combination.

### `validate.job`
- Creates output directory structure in tmpdir and `$validation/$prefix/`.
- Copies R scripts, signature annotation files (signature IDs, ATC classification, compound IDs, MoA info), and IID/phenotype files to tmpdir.
- Unpacks all `.drugs.tar.gz` files from `$drugs_pred/$prefix/` into a local `drugs/` directory.
- If `enrich=T`: runs `validate.enrich.R`; tars and copies enrichment results.
- If `ROC!=F`: runs `validate.roc.R`.

### `validate.enrich.R`
- For a given ATC code, evaluates whether drug predictions are enriched for that class.
- Loads all per-individual `.drugs` files for the cases; subsets signatures to those overlapping the ATC database.
- **Individual-level enrichment:** for each ATC sub-class, tests (Fisher's exact, one-sided) whether significantly predicted drugs (p < 0.05) are over-represented among that class's drugs. Saves per-individual OR and p-value.
- **Joint enrichment:** pools significant/non-significant hits across all cases, runs Fisher test per sub-class. Outputs results ordered by p-value, annotated with class names, drug lists, and mechanism of action.

### `validate.roc.R`
- Computes per-individual AUC values to assess drug prediction accuracy.
- "Pure cases" = individuals whose medication belongs to exactly one ATC class at the given level (no overlap with other classes).
- Optionally subsets drug signatures to specific cell types (CNS + normal blood).
- Evaluates predictions at multiple MoA granularity levels: `moa → broad → cat → cat_mono → cat_mono_exin`.
- For each individual: reads their `.drugs` file, identifies their actual drugs (from ATC class), sweeps p-value thresholds to compute TPR/FPR, calculates AUC.
- Saves per-individual AUC files and median AUC summary tables.

---

## File locations (defined in `settings.sh`)

| Variable | Contents |
|---|---|
| `$grex/all_ind/$tissue/` | Per-gene per-individual GREx files from Step 1 |
| `$collated/` | Collated full and classN GREx matrices from Step 2 |
| `$drugs_pred/$prefix/` | Per-individual drug prediction results from Step 3 |
| `$validation/$prefix/` | ROC and enrichment results from Step 4 |
| `$sigs/` | Signature annotation files (LINCS2 IDs, ATC, MoA, compound IDs) |
| `$ukb_phenos/` | UKB phenotype files (classN IIDs, medication data) |
| `$gtex/$eqtls/$tissue/` | GTEx eQTL data used for imputation |
