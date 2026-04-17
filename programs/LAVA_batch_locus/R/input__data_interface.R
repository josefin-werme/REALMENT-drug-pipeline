#' @include global.R
#' @include input__ld_reference.R
#' @include input__sumstats.R
#' @include input__annotation.R
#' @include processing__ld_blocks.R
#' @include processing__locus.R


## composite of LD reference and any number of summary statistics inputs
## summary statistics are loaded from file according to input specification, and are aligned to the LD reference
## phenotypes are loaded on demand, or when data is explicitly accessed; phenotype info objects for all are stored in 'available'
## sampling correlations are assumed to be zero unless explicitly specified using set.sampling.correlation function
##   if partial sampling correlation matrix is set, correlations with excluded phenotypes are assumed to be zero
##   input matrix needs to have phenotype names set as column names
## data can be set to load only a subset of chromosomes
##   if trim.data=F, this only affects data loading for phenotypes with input files split by chromosome, whole-genome input files are loaded in their entirety (but other chromosomes are masked)
##   if trim.data=T, whole-genome input files are truncated to only the requested chromosomes, and the entire file will be reloaded if different chromosomes are subsequently requested
DataSet = R6::R6Class("DataSet",
	private = list(
		settings = DEFAULT.SETTINGS("DataSet", IMPORT.DEFAULTS("SummaryStatistics")),
		ld.reference = NULL,   ## LDreference object
		sum.stats = list(),   ## named list of SummaryStatistics objects for loaded phenotypes, by phenotype name
		sampling.correlation = NULL,   ## CovarianceMatrix object; if NULL, assume independence
		available = list(),   ## named list of PhenotypeInfo objects for available phenotypes, by phenotype name
		chromosomes = NULL, ## ChromosomeList object with currently requested chromosomes

		## check validity of phenotypes, match partial phenotype names against available phenotypes, and remove duplicates/expand wildcards; NB: NULL input returns NULL output
		check.phenotypes = function(phenotypes, force.load=F) {
			if (length(phenotypes) > 0) {
				phenotypes = check.phenotypes(phenotypes, names(private$available), label="input data")
				if (force.load) self$load(phenotypes)
				return(phenotypes)
			} else return (if (is.null(phenotypes)) NULL else character(0))
		},

		## load data for given phenotype and currently requested chromosomes; will do nothing if data is already loaded
		load.phenotype = function(phenotype) {
			if (!phenotype %in% names(private$available)) fatal.error("invalid phenotype provided for private function load.phenotype()")
			if (phenotype %in% names(private$sum.stats)) {
				if (!private$ld.reference$equals(private$sum.stats[[phenotype]]$get.reference())) private$sum.stats[[phenotype]] = NULL
				else if (private$sum.stats[[phenotype]]$get.chromosomes("active")$equals(private$chromosomes)) return()
			}

			if (phenotype %in% names(private$sum.stats)) {
				sum.stats = private$sum.stats[[phenotype]]; private$sum.stats[[phenotype]] = NULL ## set to NULL to maintain data integrity in case of error
				private$sum.stats[[phenotype]] = sum.stats$load(private$chromosomes)
			}	else private$sum.stats[[phenotype]] = SummaryStatistics$new(private$available[[phenotype]], private$ld.reference, private$settings)$load(private$chromosomes)
		}
	),
	public = list(
		initialize = function(ld.prefix, ..., chromosomes="all") {
			private$settings = initialize.settings(self, ...)
			private$ld.reference = plink.interface(ld.prefix, chromosomes)
			private$chromosomes = private$ld.reference$get.chromosomes("loaded")$intersection(chromosomes)
		},

		print = function(...) {
			printer = ObjectPrinter$new("DataSet")
			printer$add.parameter("LD reference", private$ld.reference$abbreviate())
			printer$add.parameter("loaded phenotypes", length(private$sum.stats))
			for (ph in names(private$sum.stats)) printer$add.line(private$available[[ph]]$abbreviate(), indent=2)

			unloaded = names(private$available)[!(names(private$available) %in% names(private$sum.stats))]
			if (length(unloaded) > 0) {
				printer$add.parameter("available phenotypes", length(unloaded))
				for (ph in unloaded) printer$add.line(private$available[[ph]]$abbreviate(), indent=2)
			}
			cat(printer$to.string())
		},

		no.pheno = function(loaded.only=T) {return(ifelse(loaded.only, length(private$sum.stats), length(private$available)))},
		phenotypes = function(loaded.only=T) {return(if (loaded.only) names(private$sum.stats) else names(private$available))},

		get.chromosomes = function() {return(private$chromosomes)},

		get.reference = function() {return(private$ld.reference)},
		get.sumstats = function(phenotypes=NULL) {
			if (!is.null(phenotypes)) {
				phenotypes = private$check.phenotypes(phenotypes, force.load=T)
				return(private$sum.stats[phenotypes])
			} else return(private$sum.stats)
		},

		## fills out to full matrix matching phenotypes argument, or full list of phenotypes in data if NULL
		get.sampling.correlation = function(phenotypes=NULL) {
			if (is.null(phenotypes)) phenotypes = names(private$sum.stats)
			else phenotypes = private$check.phenotypes(phenotypes)

			if (length(phenotypes) > 0 && !is.null(private$sampling.correlation)) return(private$sampling.correlation$expand(phenotypes, covar.value=0, variances=1, ignore.existing=T)$subset(phenotypes))
			else return(create.covariance(phenotypes))
		},

		## generate DataInterface object to provide locus-level data access; option ... argument contains settings passed on to DataInterface
		get.interface = function(phenotypes=NULL, ...) {
			if (is.null(phenotypes)) {
				if (length(private$sum.stats) == 0) input.error("no phenotypes have been loaded")
				phenotypes = names(private$sum.stats)
			} else phenotypes = private$check.phenotypes(phenotypes, force.load=T)
			if (length(phenotypes) == 0) input.error("list of phenotypes is empty")
			return(DataInterface$new(self, phenotypes, ...))
		},

		add.phenotype = function(pheno.info) {
			if (pheno.info$name() %in% names(private$available)) input.error("duplicate phenotype label '", pheno.info$name(), "'; please rename")

			localized = unique(unlist(lapply(c(pheno.info, private$available), function(ii) {ii$get.meta("localized")$names})))
			multi = unique(unlist(lapply(c(pheno.info, private$available), function(ii) {ii$get.meta("multi")$name})))
			IF.ANY(intersection = intersect(multi, localized), THEN=input.error("cannot add phenotype '", pheno.info$name(), "', %parameter% ", items.and=intersection, " %occur[-s]% as both localized phenotype unit and multi-input type"))

			private$available[[pheno.info$name()]] = pheno.info
			invisible(self)
		},

		## chromosomes argument is one or a vector of numeric chromosome codes (chromosome X can be denoted as 'X' or 23), or 'all'
		## forces immediate reload of already loaded phenotypes, unless no change to reference data and phenotype has all required chromosomes already loaded
		set.chromosomes = function(chromosomes="all") {
			chromosomes = validate.chromosomes(chromosomes)
			if (!chromosomes$equals(private$chromosomes)) {
				available = private$ld.reference$get.chromosomes("available"); missing = chromosomes$difference(available)
				if (!missing$is.empty() && !(chromosomes$has.all() && missing$equals("X"))) throw.warning("%chromosome% ", missing$to.string(), " %is% not available in genotype reference data")

				required = chromosomes$intersection(available)
				if (!private$ld.reference$get.chromosomes("loaded")$contains(required)) private$ld.reference = private$ld.reference$set.chromosomes(required)
				private$chromosomes = required

				for (ph in names(private$sum.stats)) private$load.phenotype(ph)
			}
			invisible(self)
		},

		set.sampling.correlation = function(covar) {
			private$sampling.correlation = check.covariance(covar, allow.NA=T, type="sampling correlation matrix")$standardize()
			invisible(self)
		},

		## load specified phenotypes, or all if none specified
		load = function(phenotypes=NULL) {
			if (is.null(phenotypes)) phenotypes = names(private$available)
			else phenotypes = private$check.phenotypes(phenotypes)

			for (ph in phenotypes) private$load.phenotype(ph)
			invisible(self)
		},

		## unload specified phenotypes, or all if none specified
		## NB: memory may not be freed immediately if object is still referenced elsewhere
		unload = function(phenotypes=NULL) {
			if (!is.null(phenotypes)) {
				phenotypes = private$check.phenotypes(phenotypes)
				private$sum.stats[phenotypes] = NULL
			} else private$sum.stats = list()
			invisible(self)
		}
	)
)



## provides an interface to retrieve LD and summary statistics for specific phenotypes and locus
## phenotypes are set on initialization, the locus and phenotype subunits are set and updated later

## NB: DataInterface is directly dependent on objects stored in DataSet, and nay become invalid if changes are made in DataSet
## DataInferface objects are meant to be transient, discard after use
## note that LocusData/ProcessedLocus objects generated by DataInterface do exist independently, and do not depend on any previously created volatile objects

## process() will emit a ProcessedLocus (or LocusData) object, containing all the summary statistics and additional data for the current locus
## - statistics are filtered down to SNPs shared by (all units of) all loaded phenotypes (except truncated ones), after removing any that do not meet the minimum.snps criterium
## - after setting the SNP selection, phenotypes are filtered a second time on minimum.snps (and minimum.components); this does NOT trigger any subsequent change in the SNP selection
## a locus must be set by supplying a LocusDefinition object via set.locus(),
## if internal annotation is set (see set.annotation() and create.annotation() functions), set.locus() also accepts locus names and numeric indices
## - in this case, the next.locus() function can also be used to iterate over loci
## if the data contains localized or multi-input phenotypes, by default all units for which sufficient SNPs are present in the locus will be included
## - use the set.units() function to select specific unit values only; note that the unit selection will be reset when changing the locus

## STATUS: volatile, deep (LDreferenceData, SummaryStatistics)
DataInterface = R6::R6Class("DataInterface",
	private = list(
		settings = DEFAULT.SETTINGS("DataInterface",
			IMPORT.DEFAULTS("LDblock", "LocusData"),
			correlation.bound = add.info(globals$sampling.correlation.truncation, type="numeric", min=0),  ## threshold on sampling correlations to round them to 0
			load.individuals = FALSE  ## if TRUE, loads IDs and other data from .fam file and propagates them to LocusData objects
		),
		ld.reference = NULL, ## LDreferenceData object
		sum.stats = NULL, ## named list of SummaryStatistics object per phenotype
		sampling.correlation = NULL, ## CovarianceMatrix object
		units = list(
			unit.index = list(), ## named list of phenotypes per unit
			pheno.index = list(), ## named list of units per phenotype
			types = list(), ## names list of type per unit (localized/multi)
			values = list() ## named list of available unit values per unit name
		),
		annotation = NULL, ## active LocusAnnotation object (if any)
		current = list(
			locus = NULL, ## LocusDefinition object
			units = list(), ## named list of specified unit values per unit name (if NULL, use all available in locus)
			data = NULL, ## list containing cached summary statistics loaded by load.data
			ld.block = NULL ## cached LDblock object
		),

		## loads summary statistics per phenotype/units, and computes counts
		load.data = function() {
			if (is.null(private$current$locus)) input.error("no locus has been set")

			if (length(private$current$data) == 0) {
				statistics = list(); counts = list()

				## load in summary statistics and get counts per unit (if any)
				snp.index = private$current$locus$map.index(private$ld.reference$get.snp.info())
				for (ph in names(private$sum.stats)) {
					is = add.names(private$sum.stats[[ph]]$is("localized", "multi"), c("localized", "multi"))
					if (any(is)) {
						multi.name = private$sum.stats[[ph]]$get.units(include="multi", add.values=F)$names
						unit.values = add.names(lapply(private$units$pheno.index[[ph]], function(u) {private$current$units[[u]]}), names=private$units$pheno.index[[ph]])
						locus.stats = private$sum.stats[[ph]]$get.statistics(unit.values, select.ids=snp.index, multi.mode="list")

						if (is["localized"]) {
							if (is["multi"]) unit.values = unit.values[names(unit.values) != multi.name]
							else locus.stats = list(locus.stats)

							stats = list(); values = NULL; unit.fill = sapply(unit.values, is.null)
							for (i in seq_along(locus.stats)) {
								curr.stats = locus.stats[[i]]
								if (any(unit.fill)) unit.values[unit.fill] = lapply(curr.stats[names(unit.values)[unit.fill]], unique)
								value.sets = combinations(unit.values)
								for (j in 1:nrow(value.sets)) stats = list.append(stats, curr.stats[column.equals(curr.stats, value.sets[j,,drop=F]),])
								if (is["multi"]) {value.sets = cbind(names(locus.stats)[i], value.sets); names(value.sets)[1] = multi.name}
								values = rbind(values, value.sets)
							}
							locus.stats = stats
						} else {values = data.frame(names(locus.stats), stringsAsFactors=F); names(values)[1] = multi.name}

						if (length(locus.stats) > 0) {
							for (i in seq_along(locus.stats)) {
								curr.stats = locus.stats[[i]]
								if (nrow(curr.stats) >= private$settings$get("minimum.snps")) {
									curr.info = private$sum.stats[[ph]]$get.info()$instantiate.units(values[i,,drop=F])
									statistics[[curr.info$name()]] = list(phenotype=ph, info=curr.info, data=curr.stats[!(names(curr.stats) %in% names(values))])
								}
							}

							values$no.snps = sapply(locus.stats, nrow)
							counts[[ph]] = values
						}
					} else {
						curr.stats = private$sum.stats[[ph]]$get.statistics(select.ids=snp.index)
						if (nrow(curr.stats) >= private$settings$get("minimum.snps")) statistics[[ph]] = list(phenotype=ph, info=private$sum.stats[[ph]]$get.info(), data=curr.stats)
					}
				}

				available = add.colnames(matrix(F, nrow=snp.index$size(), ncol=length(statistics)), names(statistics))
				for (i in seq_along(statistics)) available[,i] = snp.index$snps() %in% statistics[[i]]$data$internal.id
				private$current$data = list(snp.index=snp.index, available=available, statistics=statistics, counts=counts)
			}
			invisible(private$current$data)
		},

		## prune down data to shared SNPs (excluding truncated phenotypes)
		filter.statistics = function(data, phenotypes=NULL) {
			if (!is.null(phenotypes)) {
				include = list.extract(data$statistics, "phenotype", flatten=T) %in% phenotypes
				if (!all(include)) {
					data$statistics = data$statistics[include]
					data$available = data$available[,include]
				}
			}

			truncated = sapply(data$statistics, function(s) {s$info$is("truncated")})
			if (any(truncated)) {
				if (!all(truncated)) {
					if (sum(!truncated) > 1) keep = apply(data$available[,!truncated], 1, all)
					else keep = data$available[,!truncated]
					data$statistics = data$statistics[apply(data$available[keep,], 2, sum) >= private$settings$get("minimum.snps")]
				} else keep = apply(data$available, 1, any)
			} else keep = apply(data$available, 1, all)
			data$snp.index = SNPindex$new(data$snp.index$snps()[keep])
			return(data[c("snp.index", "statistics")])
		}
	),
	public = list(
		initialize = function(data, phenotypes, ...) {
			private$settings = initialize.settings(self, ...)
			private$ld.reference = data$get.reference()
			private$sum.stats = data$get.sumstats(phenotypes)
			private$sampling.correlation = data$get.sampling.correlation(phenotypes)

			if (any(duplicated(names(private$sum.stats)))) input.error("duplicate phenotypes selected in DataInterface")

			has.units = sapply(private$sum.stats, function(ss) {any(ss$is("localized", "multi"))})
			if (any(has.units)) {
				for (ph in names(private$sum.stats[has.units])) {
					curr = private$sum.stats[[ph]]$get.units(add.values=T)
					private$units$pheno.index[[ph]] = curr$names
					private$units$types = list.merge(private$units$types, curr$types)
					for (u in curr$names) {
						private$units$unit.index[[u]] = c(private$units$unit.index[[u]], ph)
						private$units$values[[u]] = unique(c(private$units$values[[u]], curr$values[[u]]))
					}
				}
			}
		},

		print = function(...) {
			printer = ObjectPrinter$new("DataInterface")
			printer$add.parameter("LD reference", private$ld.reference$abbreviate())
			printer$add.list("no. phenotypes", lapply(private$sum.stats, function(ss) {ss$get.info()$abbreviate()}), values.only=T)

			if (!is.null(private$current$locus)) printer$add.line("current locus:", private$current$locus$name())
			else printer$add.line("current locus: NA")
			if (length(names(private$units$unit.index)) > 0) {
				printer$add.line("current units:")
			 	for (unit in names(private$units$unit.index)) {
			 		values = ifelse(is.null(private$current$units[[unit]]), "[ALL]", paste(private$current$units[[unit]], collapse=", "))
			 		printer$add.parameter(paste0(unit, " (", private$units$types[[unit]], ")"), values, indent=2)
			 	}
			}
			printer$add.list("settings", private$settings$get.all(), show.size=F)

			cat(printer$to.string())
		},

		## installs internal annotation object
		## allows set.locus() to take locus name and numeric index arguments, and iteration over all loci using next.locus()
		set.annotation = function(annotation=NULL) {
			if (!is.null(annotation)) {
				check.types(annotation="LocusAnnotation")
				self$set.locus(NULL)
			}
			private$annotation = annotation
			invisible(self)
		},

		## create annotation based on a localized phenotype unit, and set as internal annotation
		## if multiple phenotypes have the same unit, annotations will be created separately for each and merged (taking either their union or intersection)
		create.annotation = function(unit.name, phenotypes=NULL, merge.mode=c("union", "intersection")) {merge.mode = match.arg(merge.mode)
			if (length(unit.name) > 1) input.error("cannot create multi-unit annotation")
			if (!(unit.name %in% names(private$units$unit.index))) input.error("unknown unit ", unit.name)
			if (private$units$types[unit.name] == "multi") input.error("cannot create annotation for multi-input type '", unit.name, "'")

			if (length(phenotypes) > 0) {
				phenotypes = check.phenotypes(phenotypes, names(private$sum.stats), label="data interface")
				matched = phenotypes %in% names(private$units$unit.index)
				if (!all(matched)) input.error("unit '", unit.name, "' not available for %phenotype% ", items.or=phenotypes[!matched])
			} else phenotypes = private$units$unit.index[[unit.name]]

			annot = private$sum.stats[[phenotypes[1]]]$create.annotation(unit.name)
			if (length(phenotypes) > 1) {
				for (p in 2:length(phenotypes)) annot = annot$merge(private$sum.stats[[phenotypes[p]]]$create.annotation(unit.name), merge.mode=merge.mode)
			}

			self$set.annotation(annot)
			invisible(self)
		},

		## advance to next locus in internal annotation, return FALSE if no further remaining
		next.locus = function() {
			if (is.null(private$annotation)) input.error("no annotation has been specified for data interface")
			if (!is.null(private$current$locus)) {
				if (private$current$locus$annot.id() != private$annotation$get.id()) {
					throw.warning("current locus does not belong to internal annotation, cannot iterate; resetting current locus")
					next.id = Inf
				} else next.id = private$current$locus$id() + 1
			} else next.id = 1

			if (next.id <= private$annotation$size()) {self$set.locus(next.id); invisible(TRUE)	}
			else {self$set.locus(NULL); invisible(FALSE)}
		},

		## provide LocusDefinition object obtained from a LocusAnnotation (output of read.loci() function)
		## if internal annotation is set, locus.info can also be string or numeric index, to retrieve corresponding locus from internal annotation
		## if detect.units=T, will automatically fill in unit values for units with name matching the LocusDefinition type
		set.locus = function(locus.info, detect.units=T) {
			if (!is.null(locus.info) && !has.type(locus.info, "LocusDefinition")) {
				if (!is.null(private$annotation)) {
					locus.info = private$annotation$get.locus(locus.info)
				} else input.error("argument 'locus.info' must have type 'LocusDefinition' unless annotation has been specified for data interface")
			}
			private$current = list(locus = locus.info, units = list())

			if (!is.null(locus.info) && detect.units && length(names(private$units$unit.index)) > 0 && !("unknown" %in% locus.info$type())) {
				values = list()
				for (type in locus.info$type()[locus.info$type() %in% names(private$units$unit.index)]) values[[type]] = locus.info$name()
				if (length(values) > 0) self$set.units(values)
			}

			invisible(self)
		},

		## specify unit values as (list of) named arguments (can be vectors of multiple values)
		## if no arguments are provided, the unit selection is reset; changing the locus will also reset the unit selection
		set.units = function(...) {
			units = flatten.arglist(...)
			if (length(units) > 0) {
				if (is.null(names(units))) input.error("unit names not specified")
				IF.ANY(unknown = !(names(units) %in% names(private$units$unit.index)), THEN=input.error("unknown %unit% ", items.and=names(units)[unknown]))
			}

			private$current$data = private$current$units = list()  ## clear cached statistics and previously set units
			for (u in names(units)) {
				values = unique(units[[u]])
				IF.ANY(invalid = !(values %in% private$units$values[[u]]), THEN=input.error("unknown %value% ", items.and=values[invalid], " for unit '", u, "'"))
				private$current$units[[u]] = values
			}
			invisible(self)
		},

		available.units = function() {
			units = names(private$units$unit.index); out = NULL
			if (length(units) > 0) {
				counts = lapply(private$load.data()$counts, function(curr) {curr[units[!(units %in% names(curr))]] = NA; return(curr)})
				for (ph in names(counts)) out = rbind(out, data.frame(phenotype=ph, counts[[ph]][,c(units, "no.snps")], stringsAsFactors=F))
			}
			return(out)
		},

		## retrieves and processes LD and summary statistics, returns LocusData object
		## if discard.empty=T, discard phenotypes with insufficient SNPs in input and return reduced data object, rather than terminating
		## if compute.marginal=T, process immediately into ProcessedLocus object, otherwise return LocusData object
		## the ... argument takes an optional list of (partial) phenotype names to include; if not specified, all phenotypes are loaded
		process = function(..., discard.empty=T, compute.marginal=T) {
			phenotypes = if (length(c(...)) > 0) check.phenotypes(c(...), names(private$sum.stats), label="data interface") else names(private$sum.stats)
			if (length(phenotypes) == 0) input.error("no phenotypes have been specified")
			data = private$filter.statistics(private$load.data(), phenotypes=phenotypes)

			## verify configuration
			empty.pheno = !(phenotypes %in% list.extract(data$statistics, "phenotype", flatten=T))
			if (all(empty.pheno) || data$snp.index$size() < private$settings$get("minimum.snps")) data.error("current configuration contains fewer than ", private$settings$get("minimum.snps"), " SNPs, cannot load data")
			if (any(empty.pheno)) {
				msg = list("locus contains summary statistics for fewer than ", private$settings$get("minimum.snps"), " SNPs for %phenotype% ", items.and=phenotypes[empty.pheno])
				if (discard.empty) {msg = c(msg, "; discarding %phenotype%"); throw.warning(msg) } else data.error(msg)
			}

			sample.sizes = unlist(lapply(data$statistics, function(s) {ifelse("sample.size" %in% names(s$data), mean(s$data$sample.size), s$info$get.global("sample.size"))}))
			if (any(is.na(sample.sizes))) fatal.error("sample size information is not available for all phenotypes")
			max.components = max(globals$local.minimum.components, globals$component.sample.ratio * min(sample.sizes))


			## compute LD block, or retrieve existing if SNP selection is unchanged
			if (is.null(private$current$ld.block) || !private$current$ld.block$get.snp.index()$equals(data$snp.index)) {
				ld.block = private$ld.reference$get.block(data$snp.index, load.individuals=private$settings$get("load.individuals"), max.components=max.components, private$settings)
				if (ld.block$no.components() < private$settings$get("minimum.components")) data.error("fewer than ", private$settings$get("minimum.components"), " genetic principal components in locus genotype data")
				private$current$ld.block = ld.block
			} else ld.block = private$current$ld.block


			## filter statistics to SNP selection, and create SumstatBlocks
			stat.blocks = list(); truncate.failed = c()
			for (ph.out in names(data$statistics)) {
				curr.stats = data$statistics[[ph.out]]$data; curr.ld = ld.block
				if (data$statistics[[ph.out]]$info$is("truncated")) {
					curr.ld = ld.block$intersection(curr.stats$internal.id)
					if (curr.ld$no.components() < private$settings$get("minimum.components")) {truncate.failed = c(truncate.failed, ph.out); next}
					curr.stats = tryCatch(curr.ld$align.data(curr.stats, mode="truncate"), error=function(err) {private$fatal("invalid input to SumstatBlock object, ", err)})
					if (nrow(curr.stats) < private$settings$get("minimum.snps")) {truncate.failed = c(truncate.failed, ph.out); next}
				} else curr.stats = tryCatch(curr.ld$align.data(curr.stats, mode="truncate"), error=function(err) {private$fatal("invalid input to SumstatBlock object, ", err)})

				BlockType = if (data$statistics[[ph.out]]$info$is("binary")) BinarySumstatBlock else ContinuousSumstatBlock
				stat.blocks[[ph.out]] = BlockType$new(curr.ld, curr.stats, data$statistics[[ph.out]]$info)
			}

			if (!is.null(truncate.failed)) {
				failed.index = names(data$statistics) %in% truncate.failed
				failed.pheno = unique(list.extract(data$statistics[failed.index], "phenotype", flatten=T))
				empty.pheno = !(failed.pheno %in% list.extract(data$statistics[!failed.index], "phenotype", flatten=T))
				if (any(empty.pheno)) {
					msg = list("not enough SNPs and/or genetic principal components in locus genotype data for truncated %phenotype% ", items.and=failed.pheno[empty.pheno])
					if (discard.empty) {msg = c(msg, "; discarding %phenotype%"); throw.warning(msg) } else data.error(msg)

				}
				data$statistics = data$statistics[!failed.index]
				if (length(data$statistics) == 0) data.error("after filtering on minimum number of SNPs and genetic principal components, no phenotypes left to analyse")
			}


			## create sampling correlation matrix
			margin = list.extract(data$statistics, "phenotype", flatten=T)
			corr.matrix = add.dimnames(private$sampling.correlation$get()[margin,margin,drop=F], names(data$statistics))
			if (any(duplicated(margin))) corr.matrix[outer(margin, margin, FUN="==") & row(corr.matrix) != col(corr.matrix)] = NA
			if (length(margin) > 1) {
				truncate = !is.na(corr.matrix) & row(corr.matrix) != col(corr.matrix) & abs(corr.matrix) < private$settings$get("correlation.bound")
				corr.matrix[truncate] = 0
			}
			corr.matrix = check.covariance(corr.matrix, allow.NA=T, type="sampling correlation matrix")


			locus.data = LocusData$new(private$current$locus, ld.block, stat.blocks, corr.matrix, private$settings)
			if (compute.marginal) return(locus.data$process())
			else return(locus.data)
		},

		phenotypes = function() {return(names(private$sum.stats))},

		get.settings = function(as.list=F) {return(if (as.list) private$settings$get.all() else private$settings$clone())},
		get.input = function() {return(list(ld.reference=private$ld.reference, sum.stats=private$sum.stats, sampling.correlation=private$sampling.correlation))},

		get.annotation = function() {return(private$annotation)},
		get.locus = function() {return(private$current$locus)},
		get.units = function() {return(list(units=private$units$unit.index, values=private$current$units))}
	)
)




