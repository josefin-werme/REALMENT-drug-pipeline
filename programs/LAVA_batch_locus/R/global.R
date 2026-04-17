## keep these in order to main correct load order
#' @include utils__logging.R
#' @include utils__input_validation.R
#' @include utils__data_input.R
#' @include utils__list.R
#' @include utils__matrix.R
#' @include utils__data_frame.R
#' @include utils__math.R
#' @include utils__misc.R

.onLoad = function(libname, pkgname) {
	packageStartupMessage(paste("Running LAVA version", packageVersion("LAVA")))
}

set.mode = function(...) {
	modes = unique(c(...)); unknown = !(modes %in% c("normal", "beta", "dev"))
	if (any(unknown)) fatal.error("unknown %mode% ", modes[unknown])

	if (length(modes) > 1) modes = modes[modes != "normal"]
	globals$mode = modes
}

## can be multiple
check.mode = function(...) {any(c(...) %in% globals$mode)}


globals = as.environment(list(
	mode = "normal",
	warning.message.state = F,

	log.indent.string = "...",   ## string representing single indent step in message log

	min.pvalue = 1e-300,   ## lower bound on input p-values (and corresponding bound on input test statistics)

	ref.minimum.N = 50,   ## minimum sample size for reference data
	sampling.correlation.truncation = 0.01,   ## minimum absolute value, below which sampling correlation values are treated as zero

	locus.minimum.snps = 50,   ## minimum number of SNPs for processing locus
	locus.minimum.components = 10,   ## minimum number of PCs for processing locus
	locus.minimum.lower.bound = 1,   ## jw: changed to 1 ## lower bound on what minimum SNPs and components can be set to by user

	pca.prune.threshold = 0.99,   ## proportion of variance to retain when pruning principal components
	component.sample.ratio = 0.75,   ## maximum proportion of PCs relative to sample size
	h2.component.ratio = 0.1,   ## cutoff for ratio of components to sample size, setting h2 estimates to NA if exceeded

	stat.direction.avg.threshold = 0.05,   ## SNP effect direction proportion cutoff for direction warning (proportion below value or above 1 - value)
	stat.deflation.threshold = 0.25,   ## SNP statistic SD cutoff for deflation warning

	variant.match.warning.proportion = 0.01,   ## proportion of SNPs in input matched to reference below which to issue warning about lack of overlap
	variant.overlap.warning.proportion = 0.05,   ## proportion of SNPs overlapping between phenotype pairs below which to issue warning about lack of overlap
	variant.overlap.failure.count = 10,   ## absolute count of SNPs overlapping between phenotype pairs below which to throw error

	rg.invalid.bounds = 1.25,   ## default local rg absolute bound (considered invalid if exceeded)
	regression.invalid.bounds = 1.5,   ## default standardized regression coefficient bound (considered invalid if exceeded)
	partial.cor.max.r2 = 0.95,   ## default maximum explained variance for univariate models in partial correlation analysis
	conf.interval.threshold = 0.05,  ## p-value threshold below which to compute confidence intervals (if computation is set to be conditional on p-values)

	decomposition.eigen.threshold = 1e-4,   ## threshold for relative eigenvalues to determine invertability
	adaptive.pval.threshold = c(1e-4, 1e-6),   ## default thresholds for adaptive p-value procedure
	analysis.signif.decimals = 6,   ## number of significant decimals to round output values to

	simulation.minimum.N = 100,   ## minimum sample size for simulations

	plural.form = list(is="are", was="were", does="do", has="have", directory="directories", property="properties"),   ## dictionary for irregular singular to plural conversion in log/error messages

	## repository of default settings stored by class name, specified in other .R source files
	## add new defaults by calling DEFAULT.SETTINGS(class.name, ...)
	##   returns NULL, can use settings = DEFAULT.SETTINGS(class.name, ...) inside private part of class definition (still needs to be initialized in constructor!)
	##   include IMPORT.DEFAULTS(class.name, ...) token with one or more class names to automatically include defaults from other classes (inserted before own defaults)
	default.settings = DefaultSettings$new(),


	## automatically recognized headers for different file types, mapping onto internal parameter names (see InputInterface class in input__core.R)
	##   column names are matched case insensitive, and are converted internally to corresponding parameter names
	##   parameter names are included as valid column names as well, even if not explicitly listed here

	## locus annotation files
	annot.header.index = list(
		locus = c("LOC", "locus", "name", "gene"),
		chromosome = c("CHR", "CHROM", "chromosome"),
		start = c("start", "from"),
		stop = c("stop", "end", "to"),
		snps = c("SNPS", "IDs")
	),

	## phenotype info files
	info.header.index = list(
		phenotype = c("phenotype", "pheno"),
		filename = c("filename", "file.name", "file"),
		sample.size = c("N", "NMISS", "N_analyzed", "sample.size"),
		no.cases = c("cases", "no.cases", "ncase", "ncas"),
		no.controls = c("controls", "no.controls", "ncontrol", "ncon"),
		case.proportion = c("case.prop", "prop.case", "case.proportion"),   ## sample level
		prevalence = c("prevalence", "prev"),   ## population level
		unit = c("unit"),   ## unit specifiers for localized phenotypes
		input.type = c("type", "data.type", "input.type"),   ## currently: binary, continuous
		properties = c("properties", "prop", "traits"),   ## miscellaneous additional properties
		parameters = c("paramaters", "params")   ## user-defined parameter=column mappings for sumstats file
	),
	phenotype.properties = c("truncated", "multi", "prealigned"),   ## recognized phenotype properties for 'properties' column

	## summary statistic input files
	sumstats.header.index = list(
		snp.id = c("SNP", "ID", "SNPID_UKB", "SNPID", "MarkerName", "RSID", "RSID_UKB", "snp.id"),
		beta = c("B", "BETA"),
		odds.ratio = c("OR"),
		log.odds = c("logOR", "logOdds"),
		statistic = c("Z", "T", "ZSTAT", "STAT", "Zscore", "statistic"),
		p.value = c("P", "Pval", "Pvalue", "p.value"),
		allele1 = c("A1", "ALT", "allele1"),
		allele2 = c("A2", "REF", "allele2"),
		sample.size = c("N", "NMISS", "N_analyzed", "sample.size"),
		no.cases = c("cases", "no.cases", "ncase", "ncas"),
		no.controls = c("controls", "no.controls", "ncontrol", "ncon"),
		case.proportion = c("case.prop", "prop.case", "case.proportion")
	)
))

DEFAULT.SETTINGS = function(class.name, ...) {globals$default.settings$add(class.name, ...)}
IMPORT.DEFAULTS = function(class.name, ...) {out = c(class.name, ...); attr(out, "IMPORT") = T; return(out)}


