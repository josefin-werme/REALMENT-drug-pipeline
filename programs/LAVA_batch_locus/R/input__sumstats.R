#' @include global.R
#' @include input__phenotypes.R
#' @include input__ld_reference.R




## reads in summary statistics and aligns them to the LD reference, filtering out variants not present in reference and entries with NA values
## validation of column specification/settings is already performed in PhenotypeInfo object, and data load is run via PhenotypeInfo object as well
## NB: if chromosomes argument is set to subset, but data stored in single input file, then full file is loaded but other chromosomes are masked

## provided counts are, computed in order
##   rows: number of variants read from input file
##   available: number remaining after filtering those not in reference
##   aligned: number remaining after filtering those that could not be aligned to reference (allele mismatch, ambiguous AT/CG)
##   valid: number remaining after filtering on missing summary statistics, other missing values in relevant columns, zero sample size, and duplicate SNP IDs
##   active: number of valid SNPs that are in active chromosome selection
## NB: if trim.data=T, first four values refer to last data load, and include specified chromosomes only; they are not updated if data is subsequently trimmed without forcing a new data load
SummaryStatistics = R6::R6Class("SummaryStatistics",
	private = list(
		settings = DEFAULT.SETTINGS("SummaryStatistics",
			snp.duplicates = c("drop", "first", "last"),   ## determines what to do with duplicate entries: discard all (drop; default), retain only the first, or only the last
			trim.data = FALSE   ## if TRUE, for files not stored by chromosome, data for chromosomes that are not needed is discarded to reduce memory and processing time (NB: this will force a reload of the whole file every time a different chromosome is requested)
		),

		pheno.info = NULL, ## PhenotypeInfo object
		ld.reference = NULL,  ## LDreference object
		sum.stats = NULL,  ## data.frame with columns: internal.id, statistic, sample.size (optional), case.proportion (binary; optional), localized unit columns (optional; stored as factor), active (optional)

		counts = list(rows=0, available=0, aligned=0, valid=0, active=0),  ## SNP counts of loaded data (see above)
		localized = list(names = NULL, counts = NULL), multi = list(name = NULL, values = NULL),  ## vector of unit names, and name/values of multi-input type

		active.snps = NULL, ## index of active SNPs
		chromosomes = list(active=NULL, loaded=NULL),  ## ChromosomeList objects (NB: this includes chromosomes for which files were missing)

		warning = function(..., filename=NULL) {throw.warning(..., .label=list(phenotype=private$pheno.info$name(), file=ifelse(is.null(filename), private$pheno.info$filename(), filename)))},
		error = function(..., filename=NULL, fatal=F) {ifelse(fatal, fatal.error, input.error)(..., .label=list(phenotype=private$pheno.info$name(), file=ifelse(is.null(filename), private$pheno.info$filename(), filename)))},

		load.data = function(chromosomes) {chromosomes = validate.chromosomes(chromosomes)
			## load in data and align
			if (private$pheno.info$by.chromosome()) {
				avail.chromosomes = chromosome.files(private$pheno.info$filename(), chromosomes, prune.missing=T)$chromosome
				if (length(avail.chromosomes) > 0) {
					for (i in seq_along(avail.chromosomes)) {
						curr = private$read.file(avail.chromosomes[i])
						if (i > 1) {
							if (identical(names(sum.stats), names(curr))) sum.stats = rbind(sum.stats, curr)
							else private$error("cannot merge chromosomes, headers are incompatible", filename=private$pheno.info$filename(avail.chromosomes[i]))
						} else sum.stats = curr
					}
				} else private$error("specified input files do not exist")
			} else sum.stats = private$read.file()

			if (private$settings$get("trim.data")) {
				sum.stats = sum.stats[abs(sum.stats$chromosome) %in% chromosomes$get(),]
				if (length(private$localized$names) > 0) sum.stats[private$localized$names] = lapply(sum.stats[private$localized$names], droplevels)
			}

			## trim SNPs not in data, or (if stored per chromosome) with chromosome code not matching file
			available = !is.na(sum.stats$internal.id)
			invalid.chr = available & (is.na(sum.stats$chromosome) | sum.stats$chromosome <= 0)
			counts = list(rows = nrow(sum.stats), available = sum(available))
			if (any(invalid.chr)) private$warning("discarding ", items.count=sum(invalid.chr), " %SNP% with mismatched chromosome %code% from input")
			sum.stats = sum.stats[available & !invalid.chr,]

			## filter out duplicates and unaligned SNPs (unaligned count includes SNPs with mismatched chromosomes)
			## NB: performing duplication check before any other filtering so selection does not depend on other statistics; chromosome mismatch filter guarantees that SNP ID cannot occur in different files
			sum.stats = private$filter.duplicates(sum.stats, dupl.columns = private$localized$names)
			if (!self$is("prealigned")) sum.stats = sum.stats[!is.na(sum.stats$direction),]
			counts$aligned = nrow(sum.stats)

			## process and validitate input values; configuration is assumed to have been validated in PhenotypeInfo object
			configuration = private$pheno.info$get.configuration()
			if ("sample.size" %in% names(sum.stats)) sum.stats$sample.size = suppressWarnings(as.numeric(sum.stats$sample.size))
			sum.stats = private$process.statistics(sum.stats, configuration$statistics)
			if (self$is("binary")) sum.stats = private$process.binary(sum.stats, configuration$metrics)
			if ("sample.size" %in% names(sum.stats)) sum.stats$sample.size[!is.finite(sum.stats$sample.size) | sum.stats$sample.size <= 0] = NA

			## discard extraneous columns (NB: chromosome column is trimmed off in set.active() later, if trim.data=T)
			multi.columns = grep("^statistic::", names(sum.stats), value=T)
			keep.columns = c("internal.id", "chromosome", "statistic", multi.columns, "sample.size", "case.proportion", private$localized$names)
			sum.stats = sum.stats[names(sum.stats) %in% keep.columns]

			is.multi = names(sum.stats) %in% multi.columns
			keep = column.none(is.na(sum.stats[!is.multi]))
			if (any(is.multi)) keep = keep & !column.all(is.na(sum.stats[!is.multi]))  ## for multi-input, only discard if all statistics are NA
			sum.stats = sum.stats[keep,]
			counts$valid = nrow(sum.stats)


			private$sum.stats = sum.stats
			private$counts = counts
			private$chromosomes$loaded = if (private$pheno.info$by.chromosome() || private$settings$get("trim.data")) chromosomes else ChromosomeList$new("all")
			private$set.active(chromosomes)
		},

		## sets specified chromosomes as active (masking/discarding, and updates relevant counts and unit values
		set.active = function(chromosomes) {
			private$active.snps = NULL
			chromosomes = validate.chromosomes(chromosomes)

			if (!private$chromosomes$loaded$equals(chromosomes)) {
				if (!("chromosome" %in% names(private$sum.stats))) private$sum.stats$chromosome = private$ld.reference$get.snp.info(private$sum.stats$internal.id, "chromosome")
				selected = private$sum.stats$chromosome %in% chromosomes$get()
				if (!all(selected)) {
					if (private$settings$get("trim.data")) {
						private$sum.stats = private$sum.stats[selected,]
						private$chromosomes$loaded = chromosomes
					} else private$active.snps = which(selected)
				}
			}

			private$sum.stats$chromosome = NULL
			private$chromosomes$active = chromosomes
			private$counts$active = ifelse(!is.null(private$active.snps), length(private$active.snps), nrow(private$sum.stats))
			if (self$is("localized")) private$localized$counts = lapply(self$get.units(add.values=T, include="localized")$values, length)
		},

		## read in data file, aligns to reference, and processes unit columns
		## - adds internal.id (NA if unmatched), chromosome (negative if inconsistent) and direction (NA if not alignable) columns (direction is not added for 'prealigned' phenotypes)
		## - processes unit columns, for localized phenotypes
		read.file = function(chromosome=NULL) {
			sum.stats = private$pheno.info$get.data(chromosome)

			if (is.null(sum.stats)) private$error("unable to load input data")
			if (nrow(sum.stats) == 0) private$error("file is empty", filename=private$pheno.info$filename(chromosome))

			## map to reference data, replacing snp.id, allele1 and allele2 columns by internal.id column (set to NA if not matched)
			## adds direction column with values -1/0/1 (set to NA for unalignable SNPs)
			## adds chromosome column (set to negative value if not matching intended chromosome)
			sum.stats = private$ld.reference$align.sumstats(sum.stats, add.chromosome=T, skip.align=self$is("prealigned"))
			if (!is.null(chromosome)) {
				mismatched = !is.na(sum.stats$chromosome) & sum.stats$chromosome != chromosome
				sum.stats$chromosome[mismatched] = -sum.stats$chromosome[mismatched]
			}

			## process columns for units, if localized
			if (self$is("localized")) {
				for (unit in private$localized$names) {
					values = sum.stats[[unit]]
					is.repeat = sum.stats[[unit]] == "."; is.repeat[1] = F
					if (any(is.repeat)) {
						index = which(!is.repeat); count = c(index[-1], length(is.repeat)+1) - index
						values = values[rep(index, times=count)]
					}
					values[values == "." | values == ""] = NA
					sum.stats[[unit]] = as.factor(values)
				}
			}
			return(sum.stats)
		},

		## select and filter out duplicate entries based on combination of internal.id and dupl.columns (if any)
		filter.duplicates = function(sum.stats, dupl.columns=NULL) {
			dupl.id = paste.columns(sum.stats, c("internal.id", dupl.columns), sep=" ")
			dupl.mode = private$settings$get("snp.duplicates")  ## drop, (keep) first, (keep) last
			is.duplicated = duplicated(dupl.id, fromLast=(dupl.mode=="last"))
			if (!is.null(dupl.columns)) is.duplicated[apply(is.na(sum.stats[dupl.columns]), 1, any)] = F
			if (any(is.duplicated)) {
				if (!(dupl.mode %in% c("first", "last"))) {
					is.duplicated = dupl.id %in% unique(dupl.id[is.duplicated])
					msg.label = "discarding all"
				} else msg.label = paste0("discarding all but ", dupl.mode)

				private$warning("found duplicate SNP IDs ", if (length(private$localized$names) > 0) "(within units)", " in input, ", msg.label, " (", items.count=sum(is.duplicated), " %line% removed)")
				sum.stats = sum.stats[!is.duplicated,]
			}
			return(sum.stats)
		},

		## creates statistic column if not yet present, and aligns to reference data (expects direction column)
		## stat.columns are assumed to be a vector of 'statistics' columns and ordered pairs of a 'p.value' column followed by a corresponding effect size column
		process.statistics = function(sum.stats, stat.columns) {
			for (stat.col in grep("^statistic", stat.columns, value=T)) sum.stats[[stat.col]] = suppressWarnings(as.numeric(sum.stats[[stat.col]]))

			## create statistic columns if not present yet
			for (i in grep("^p.value", stat.columns)) {
				p.col = stat.columns[i]
				pvals = suppressWarnings(as.numeric(sum.stats[[p.col]]))
				pvals[!is.na(pvals) & pvals < globals$min.pvalue] = globals$min.pvalue

				effect.col = stat.columns[i+1]
				if (grepl("^odds.ratio", effect.col)) {
					if (any(!is.na(sum.stats[[effect.col]]) & sum.stats[[effect.col]] < 0)) private$warning("input contains negative odds ratio values; please check that you did not provide log odds or betas")
					direction = sign(sum.stats[[effect.col]] - 1)
				}	else direction = sign(sum.stats[[effect.col]])

				avg.direction = ifelse(length(direction) > 0, mean(direction == 1), 0.5)
				if (avg.direction < globals$stat.direction.avg.threshold || avg.direction > (1-globals$stat.direction.avg.threshold)) private$warning("almost all ", gsub("::.*$", "", effect.col), " values have the same direction; please verify input is correct")

				stat.col = gsub("^p.value", "statistic", p.col)
				sum.stats[[stat.col]] = -qnorm(pvals/2) * direction
			}
			sum.stats[grep("^statistic", stat.columns, invert = T, value=T)] = NULL

			## align effect direction and validate
			for (stat.col in grep("^statistic", names(sum.stats), value=T)) {
				if (!self$is("prealigned")) sum.stats[[stat.col]] = sum.stats[[stat.col]] * sum.stats$direction

				## convert infinite-valued test statistics
				is.infinite = which(abs(sum.stats[[stat.col]]) > abs(qnorm(globals$min.pvalue/2)))
				sum.stats[[stat.col]][is.infinite] = abs(qnorm(globals$min.pvalue/2)) * sign(sum.stats[[stat.col]][is.infinite])

				if (sum(!is.na(sum.stats[[stat.col]])) > 0) {
					sq.stat.avg = mean(sum.stats[[stat.col]]^2, na.rm=T)
					if (sqrt(sq.stat.avg) < globals$stat.deflation.threshold) {
						label = if (stat.col != "statistic") list(" for '", gsub("^statistic::", "", stat.col), "'")
						private$warning("SNP summary statistics", label, " appear severely deflated; please check that input is correct")
					}
				}
			}

			return(sum.stats)
		},

		## if needed, compute sample.size and/or case.proportion columns from other input for binary phenotype
		## is expecting either just a single sample.size column, or two columns: sample.size + case.proportion, or two of sample.size, no.cases and no.controls
		process.binary = function(sum.stats, binary.columns) {
			if (length(binary.columns) >= 2) {  ## if not, either has just a sample.size column, or global values only
				for (col in binary.columns[binary.columns != "sample.size"]) sum.stats[[col]] = suppressWarnings(as.numeric(sum.stats[[col]]))
				if (!all(c("sample.size", "case.proportion") %in% binary.columns)) {
					if (!("sample.size" %in% binary.columns)) sum.stats$sample.size = sum.stats$no.cases + sum.stats$no.controls
					if ("no.cases" %in% binary.columns) sum.stats$case.proportion = sum.stats$no.cases / sum.stats$sample.size
					else sum.stats$case.proportion = 1 - sum.stats$no.controls / sum.stats$sample.size
				}
			}

			return(sum.stats)
		}
	),
	public = list(
		initialize = function(pheno.info, ld.reference, ...) {
			private$settings = initialize.settings(self, ...)
			private$pheno.info = pheno.info; private$ld.reference = ld.reference

			private$chromosomes$active = private$chromosomes$loaded = ChromosomeList$new("none")
			if (self$is("multi")) {
				multi = private$pheno.info$get.meta("multi")
				private$multi = list(name = multi$name, values = names(multi$columns))
			}
			if (self$is("localized")) private$localized$names = private$pheno.info$get.configuration()$localized
		},

		print = function(...) {
			printer = ObjectPrinter$new(class(self)[1])$parameter.list(phenotype=private$pheno.info$name(), file=private$pheno.info$filename(), type=private$pheno.info$get.traits())
			printer$add.line("variants:", private$counts$rows, "in file,", private$counts$aligned, "valid and aligned to reference")
			if (self$is("multi")) printer$add.line("multi-input: ", length(private$multi$values), " (", private$multi$name, ")", add.space=F)
			if (self$is("localized")) {for (unit in private$localized$names) printer$add.line("localized units: ", private$localized$counts[[unit]], " (", unit, ")", add.space=F)}
			cat(printer$to.string())
		},

		phenotype = function() {return(private$pheno.info$name())},
		is = function(...) {return(private$pheno.info$is(...))},  ## query phenotype traits, eg. is("binary")

		get.info = function() {return(private$pheno.info)},

		get.chromosomes = function(mode=c("active", "loaded")) {return(private$chromosomes[[match.arg(mode)]])},
		get.counts = function() {return(private$counts)},
		get.units = function(add.values=T, include=c("all", "multi","localized")) {include = match.arg(include)
			include = add.names(sapply(c("multi","localized"), function(mode) {(include == "all" || include == mode) && self$is(mode)}), c("multi","localized"))
			names = c(if (include["multi"]) private$multi$name, if (include["localized"]) private$localized$names)
			types = named.list(names=names, ifelse(names %in% private$localized$names, "localized", "multi"))

			if (add.values) {
				values = list()
				if (include["multi"]) values[[private$multi$name]] = private$multi$values
				if (include["localized"]) {
					curr.data = if (!is.null(private$active.snps)) private$sum.stats[private$localized$names][private$active.snps,] else private$sum.stats[private$localized$names]
					values = c(values, lapply(curr.data, function(f) {as.character(unique(f))}))
				}
				return(list(names=names, types=types, values=values))
			} else return(list(names=names, types=types))
		},

		get.reference = function() {return(private$ld.reference)},

		## retrieve internal IDs and unit columns (if any); convert to external IDs and additional info if add.info=T
		get.ids = function(add.info=F) {
			sum.stats = if (!is.null(private$active.snps)) private$sum.stats[private$active.snps,] else private$sum.stats
			sum.stats = sum.stats[c("internal.id", private$localized$names)]
			if (add.info) sum.stats = cbind(private$ld.reference$get.snp.info(sum.stats$internal.id, c("snp.id", "chromosome", "position", "allele1", "allele2")), subset(sum.stats, select=-internal.id))
			return(sum.stats)
		},

		## retrieve sumstats for specified internal IDs and unit values
		## ... argument specifies unit.name=values pairs (can be NULL, values can be vector of multiple)
		## the multi.mode argument determines output format for multi-input phenotypes
		## - column: output as stored with multiple columns, with 'statistics::[TYPE]' column names; still may contain NA values
		## - row: transpose into rows with additional unit column with name equal to private$multi$name and single 'statistic' column; NA values are removed
		## - list: output named list, with separate data.frame per [TYPE], and single 'statistic' column; NA values are removed
		get.statistics = function(..., select.ids=NULL, multi.mode=c("column", "row", "list")) {multi.mode = match.arg(multi.mode)
			sum.stats = if (!is.null(private$active.snps)) private$sum.stats[private$active.snps,] else private$sum.stats

			unit.values = flatten.arglist(...)
			if (length(unit.values) > 0) {
				if (!any(self$is("multi", "localized"))) private$error("unit specification has been provided for phenotype without 'localized' or 'multi' properties")

				if (self$is("multi") && private$multi$name %in% names(unit.values)) {
					multi.values = unit.values[[private$multi$name]]; unit.values = unit.values[names(unit.values) != private$multi$name]
					IF.ANY(invalid = !(multi.values %in% private$multi$values), THEN = private$error("invalid type specification for multi-input phenotype with type '", private$multi$name, "', unknown %value% ", items.and=multi.values[invalid]))
					sum.stats = sum.stats[!grepl("^statistic", names(sum.stats)) | names(sum.stats) %in% paste0("statistic::", multi.values)]
				}

				if (length(unit.values) > 0) {
					if (self$is("localized")) {
						IF.ANY(unmatched = !(names(unit.values) %in% private$localized$names), THEN = private$error("invalid unit specification for localized phenotype, unknown %unit% ", items.and=names(unit.values)[unmatched]))
						sum.stats = sum.stats[column.equals(sum.stats, unit.values),]
					} else private$error("unknown multi-input %type% ", items.and=names(unit.values))
				}
			}

			if (!is.null(select.ids)) {
				if (has.type(select.ids, "SNPindex")) select.ids = select.ids$snps()
				sum.stats = sum.stats[sum.stats$internal.id %in% select.ids,,drop=F]
			}

			if (self$is("multi") && multi.mode != "column") {
				stat.columns = grep("^statistic::", names(sum.stats), value=T); labels = as.factor(gsub("^statistic::", "", stat.columns))
				stats.base = sum.stats[!(names(sum.stats) %in% stat.columns)]

				if (multi.mode == "row") {
					stats.base[[private$multi$name]] = labels[1]; stats.base$statistic = sum.stats[[stat.columns[1]]]; stats.out = stats.base
					if (length(stat.columns) > 1) {
						for (i in 2:length(stat.columns)) {
							stats.base[[private$multi$name]] = labels[i]; stats.base$statistic = sum.stats[[stat.columns[i]]]
							stats.out = rbind(stats.out, stats.base)
						}
					}
					return(stats.out[!is.na(stats.out$statistic),])
				} else {
					stats.out = list()
					for (i in seq_along(stat.columns)) {
						curr.stats = stats.base; curr.stats$statistic = sum.stats[[stat.columns[i]]]
						stats.out[[as.character(labels[i])]] = curr.stats[!is.na(curr.stats$statistic),]
					}
					return(stats.out)
				}
			} else return(sum.stats)
		},

		## loads specified chromosomes, masking/unmasking already loaded data as needed
		## currently, for input stored by chromosome, if not all requested chromosomes are loaded, will load all requested (even if part already loaded)
		load = function(chromosomes, quiet=F) {chromosomes = validate.chromosomes(chromosomes)
			if (!private$chromosomes$loaded$contains(chromosomes)) {  ## need to load additional data
				if (!quiet) {
					log.message("reading summary statistics for phenotype '", self$phenotype(), "' (", toupper(private$pheno.info$get.traits()), ") from file ",  private$pheno.info$filename())
					if (!chromosomes$has.all() && (private$pheno.info$by.chromosome() || private$settings$get("trim.data"))) log.message("loading input for %chromosome% ", chromosomes$to.string(), .indent=2)

					config = private$pheno.info$get.configuration(add.columns=T); multi.message = NULL
					if (self$is("multi")) {
						multi.info = private$pheno.info$get.meta("multi")
						is.multi = config$parameter %in% unlist(multi.info$columns)
						multi.value = gsub("^.*::(.*)$", "\\1", config$parameter[is.multi])
						config$column[is.multi] = sapply(seq_along(multi.value), function(i) {gsub(paste0(multi.value[i], "$"), "\\[NAME]", config$column[is.multi][i], ignore.case=T)})
						config$parameter[is.multi] = gsub("^(.*)::.*", "\\1", config$parameter[is.multi])
						counts = table(config$column[is.multi])[config$column[is.multi]]
						config$parameter[is.multi] = paste0(config$parameter[is.multi], ", x", counts)
						config = unique(config)

						multi.message = list("detected summary statistics columns for ", items.count=length(multi.info$columns), " %phenotype%")
						if (multi.info$name != "multi") multi.message = c(multi.message, " of type '", multi.info$name, "'")
						if (length(multi.info$columns) <= 10) multi.message = c(multi.message, ": ", paste(names(multi.info$columns), collapse=", "))
					}
					log.message("using columns ", paste0(config$column, " (", config$parameter, ")", collapse=", "), .indent=2)
					if (!is.null(multi.message)) log.message(multi.message, .indent=2)
				}

				if (private$pheno.info$by.chromosome()) {
					missing = chromosomes$difference(private$pheno.info$get.chromosomes())
					if (!missing$is.empty()) {
						msg = list("missing input %file% for %chromosome% ", missing$to.string())
						if (missing$equals(chromosomes)) private$error(msg)
						if (!chromosomes$has.all() || missing$size(auto.only=T) > 0) private$warning(msg)
					}
				}
				private$load.data(chromosomes)

				prop.thresh = globals$variant.match.warning.proportion
				if (private$counts$available / private$counts$rows < prop.thresh) {
					if (private$counts$available == 0) private$error("no variants could be matched to reference data; please make sure the SNP ID format matches the reference data")
					else private$warning("less than ", 100*prop.thresh, "% of variants matched to reference data; please make sure the SNP ID format matches the reference data")
				}
				if (private$counts$aligned / private$counts$available < prop.thresh) {
					if (private$counts$aligned == 0) private$error("none of the variants matched to reference data could be aligned; please make sure the allele columns are correct")
					else private$warning("less than ", 100*prop.thresh, "% of variants matched to reference data could be aligned; please make sure the allele columns are correct")
				}
				if (private$counts$valid / private$counts$aligned < prop.thresh) {
					if (private$counts$aligned == 0) private$error("none of the aligned variants have valid data; please make sure your input files are correctly formatted")
					else private$warning("less than ", 100*prop.thresh, "% of aligned variants have valid data; please make sure your input files are correctly formatted")
				}
				if (!quiet) log.message("read ", private$counts$rows, " variants, of which ", private$counts$valid, " valid and aligned to reference", .indent=2)

				if (private$counts$active < private$counts$valid) {
					if (private$counts$active == 0) input.error("input data contains no aligned variants with valid data for specified %chromosome% ", chromosomes$to.string())
					if (!quiet) log.message("setting active %chromosome% ", chromosomes$to.string(), ", ", pluralize(items.count=private$counts$active, " %variant% remaining"), .indent=2)
				}
			} else if (!private$chromosomes$active$equals(chromosomes)) {  ## only need to change the active chromosomes
				if (!quiet) log.message("changing active %chromosome% to ", chromosomes$to.string(), " for phenotype '", self$phenotype(), "' (", toupper(private$pheno.info$get.traits()), ")")
				private$set.active(chromosomes)

				if (private$counts$active == 0) private$error("input data contains no aligned variants with valid data for specified %chromosome% ", .plural=chromosomes$size())
				if (!quiet) log.message("for specified %chromosome%", .plural=chromosomes$size(), ", ", pluralize(items.count=private$counts$active, " aligned %variant% with valid data are available"), .indent=2)
			} else invisible(self)

			if (!quiet && self$is("localized")) log.message("localized unit counts: ", paste0(unlist(private$localized$counts), " of type '", names(private$localized$counts), "'", collapse=", "), .indent=2)
			invisible(self)
		},

		## create locus annotation object based on unit
		create.annotation = function(unit.name) {
			if (!self$is("localized")) private$error("cannot create annotation, phenotype is not localized")
			if (length(unit.name) > 1) private$error("cannot create multi-unit annotation")
			if (self$is("multi") && unit.name == private$multi$unit) private$error("cannot create annotation for multi-input type '", unit.name, "'")
			if (!(unit.name %in% private$localized$names)) private$error("invalid unit specification for localized phenotype, unknown unit '", unit.name, "'")

			sum.stats = self$get.ids(add.info=T)
			locus.names = as.character(unique(sum.stats[[unit.name]]))
			chr = aggregate(sum.stats$chromosome, list(unit=sum.stats[[unit.name]]), function(x) {unique(x)})
			if (any(sapply(chr$x, length) != 1)) private$error("cannot create annotation for unit ", unit.name, ", some instances span multiple chromosomes")
			ranges = aggregate(sum.stats$position, list(unit=sum.stats[[unit.name]]), range)

			annot = add.names(data.frame(locus.names, chr[match(locus.names, chr$unit),2], ranges[match(locus.names, ranges$unit),2], stringsAsFactors=F), c("locus", "chromosome", "start", "stop"))
			return(RangeAnnotation$new(annot, source=private$pheno.info$filename(), type=unit.name))
		}
	)
)




