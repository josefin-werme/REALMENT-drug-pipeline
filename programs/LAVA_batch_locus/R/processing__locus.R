#' @include global.R
#' @include input__phenotypes.R
#' @include input__annotation.R
#' @include processing__ld_blocks.R
#' @include processing__univariate.R
#' @include processing__overlap.R



## combined input data for specific locus
## contains locus-level LD block, sampling correlation matrix, and list of summary statistics objects for each phenotype
## call process() to obtain corresponding ProcessedLocus object

## STATUS: stable, deep (LDblock, SumstatBlock)
LocusData = R6::R6Class("LocusData",
	private = list(
		settings = DEFAULT.SETTINGS("LocusData",
			minimum.snps = add.info(globals$locus.minimum.snps, min=globals$locus.minimum.lower.bound),
			minimum.components = add.info(globals$locus.minimum.components, min=globals$locus.minimum.lower.bound)
		),
		locus = NULL, ## LocusDefinition object
		ld.block = NULL, ## LDblock object
		sum.stats = list(), ## named list of SumstatBlock object per phenotype
		sampling.correlation = NULL, ## CovarianceMatrix object
		harmonized = list(), ## named list indicating whether sumstats for a given phenotype are aligned to the global ld.block stored above

		subset.internal = function(phenotypes) {
			private$sum.stats = private$sum.stats[phenotypes]
			private$sampling.correlation = private$sampling.correlation$subset(phenotypes)
			private$harmonized = private$harmonized[phenotypes]
			return(self)
		}
	),
	public = list(
		initialize = function(locus, ld.block, sum.stats, sampling.covariance, ...) {
			private$settings = initialize.settings(self, ...)
			if (length(sum.stats) == 0) input.error("no phenotypes provided to LocusData object")
			names(sum.stats) = sapply(sum.stats, function(ss) {ss$phenotype()})
			if (any(duplicated(names(sum.stats)))) input.error("input to LocusData object contains duplicate phenotypes")
			private$locus = locus; private$ld.block = ld.block; private$sum.stats = sum.stats

			failed = c()
			for (ph in names(private$sum.stats)) {
				private$harmonized[[ph]] = private$ld.block$equals(private$sum.stats[[ph]]$get.ld())
				if (!private$harmonized[[ph]] && !private$ld.block$get.snp.index()$contains(private$sum.stats[[ph]]$internal.id())) failed = c(failed, ph)
			}
			if (length(failed) > 0) input.error("invalid input to LocusData object, summary statistics for %phenotype% ", items.and=failed, " are not fully contained in provided LD block")

			private$sampling.correlation = check.covariance(sampling.covariance, allow.NA=T, type="sampling correlation matrix")$standardize()
			if (private$sampling.correlation$size() != length(sum.stats) || !all(private$sampling.correlation$get.names() == names(sum.stats))) input.error("sampling correlation matrix input to LocusData object does not match input summary statistics")
		},

		print = function(...) {
			printer = ObjectPrinter$new("LocusData")
			printer$add.parameter("locus", private$locus$name())$add.parameter("no. snps", self$no.snps())$add.parameter("no. components", self$no.components())
			printer$add.list("no. phenotypes", lapply(private$sum.stats, function(ss) {ss$phenotype.info()$abbreviate()}), values.only=T)
			cat(printer$to.string())
		},

		process = function() {return(ProcessedLocus$new(self, private$settings))},

		phenotypes = function() {return(names(private$sum.stats))},

		no.pheno = function() {return(length(private$sum.stats))},
		no.snps = function(total=T) {return(if (!total) unlist(lapply(private$sum.stats, function(s) {s$no.snps()})) else private$ld.block$no.snps())},
		no.components = function(total=T) {return(if (!total) unlist(lapply(private$sum.stats, function(s) {s$no.components()})) else private$ld.block$no.components())},

		get.locus = function() {return(private$locus)},
		get.ld = function() {return(private$ld.block)},
		get.sumstats = function() {return(private$sum.stats)},
		get.correlations = function() {return(private$sampling.correlation)},

		is.harmonized = function(aggregate=T) {return(if (!aggregate) unlist(private$harmonized) else all(unlist(private$harmonized)))},

		## keep only specified phenotypes if include=T, discard specified phenotypes otherwise
		subset = function(phenotypes, include=T, empty.mode=c("warning", "error", "quiet")) {empty.mode = match.arg(empty.mode)
			if (length(phenotypes) == 0) phenotypes = character(0)
			phenotypes = check.phenotypes(phenotypes, self$phenotypes(), label="locus data")
			if (!include) phenotypes = self$phenotypes()[!(self$phenotypes() %in% phenotypes)]

			if (length(phenotypes) == 0 && self$no.pheno() > 0 && empty.mode != "quiet") ifelse(empty.mode == "warning", throw.warning, data.error)("no phenotypes remaining after subsetting")
			return(self$clone()$.__enclos_env__$private$subset.internal(phenotypes))
		}
	)
)



## modifiable version for use with simulations
## can update SumstatBlock objects and instance IDs of phenotypes

## STATUS: volatile, deep (LDblock, SumstatBlock)
SimulationLocusData = R6::R6Class("SimulationLocusData",
	inherit = LocusData,
	private = list(),
	public = list(
		initialize = function(locus, ld.block, phenotypes, sampling.covariance, ...) {
			if (length(phenotypes) == 0 || !all(sapply(phenotypes, function(pi) {has.type(pi, "SimulatedPhenotype")}))) input.error("invalid phenotype input to SimulationLocusData object")
			sampling.covariance = sampling.covariance$subset(sapply(phenotypes, function(pi) {pi$name()}))  ##align matrix to order of phenotypes, for later updates
			sum.stats = list()
			stats.tpl = data.frame(internal.id = ld.block$get.snp.index()$snps(), statistic=NA)

			for (i in seq_along(phenotypes)) {
				BlockType = if (phenotypes[[i]]$is("binary")) BinarySumstatBlock else ContinuousSumstatBlock
				sum.stats[[i]] = BlockType$new(ld.block, stats.tpl, phenotypes[[i]])
			}

			super$initialize(locus, ld.block, sum.stats, sampling.covariance, ...)
		},

		update.id = function(id) {
			for (i in seq_along(private$sum.stats)) private$sum.stats[[i]]$.__enclos_env__$private$pheno.info = private$sum.stats[[i]]$phenotype.info()$set.id(id)
			names(private$sum.stats) = names(private$harmonized) = sapply(private$sum.stats, function(ss) {ss$phenotype()})
			private$sampling.correlation = private$sampling.correlation$set.names(names(private$sum.stats))
			invisible(self)
		},

		set.statistics = function(stats) {
			check.types(stats="data.frame")
			data.names = cbind(names(private$sum.stats), sapply(private$sum.stats, function(ss) {ss$phenotype.info()$name(base.only=T)}))
			index = sapply(names(stats), function(ph) {out = which(data.names[,1] == ph | data.names[,2] == ph); return(ifelse(length(out) == 1, out, NA))})
			IF.ANY(unknown = is.na(index), THEN=input.error("updated statistics provided for unknown %phenotype% ", items.and=names(stats)[unknown]))
			IF.ANY(missing = which(!(1:self$no.pheno() %in% index)), THEN=input.error("no updated statistics provided for %phenotype% ", names(private$sum.stats)[missing]))
			for (i in seq_along(index)) {
				curr.stats = private$sum.stats[[index[i]]]; curr.stats$clear()
				curr.stats$.__enclos_env__$private$snp.data$statistic = stats[,i]
			}
			invisible(self)
		}
	)
)


## combined marginal statistics for a locus
## contains locus-level LD block, sampling correlations, and list of processed summary statistics per phenotype
## processing status can be: available, negative, failed
## - failed: a processing error occurred, or the sigma estimate is negative or zero
## - negative: processing succeeded, but h2 / omega variance estimate is negative or zero
## - positive: same as 'available'
## - processed: inverse of 'failed'
## summarize() function provides an overview of the status of all phenotypes, as well as marginal statistics (univariate p-value, h2 estimates)

## STATUS: stable, deep (LDblock, MarginalEstimates)
ProcessedLocus = R6::R6Class("ProcessedLocus",
	private = list(
		settings = DEFAULT.SETTINGS("ProcessedLocus", IMPORT.DEFAULTS("LocusData")),
		phenotype.info = list(
			info = list(),	## named list of PhenotypeInfo objects
			status = list() ## named list with status (failed/negative/available) per phenotype
		),
		data = list(
			locus = NULL,  ## LocusDefinition object
			ld.block = NULL,  ## LDblock object
			correlations = NULL  ## sampling correlation matrix (CovarianceMatrix object), subset to non-failed phenotypes
		),
		estimates = list(), ## named list of MarginalEstimates object per phenotype (excluding failed)

		warning = function(..., phenotype=NULL) {throw.warning(..., .label=list(phenotype=phenotype, "locus ID"=private$data$locus$name()))},
		failure = function(..., phenotype=NULL) {data.error(..., .label=list(phenotype=phenotype, "locus ID"=private$data$locus$name()))},
		error = function(..., phenotype=NULL, fatal=F) {ifelse(fatal, fatal.error, input.error)(..., .label=list(phenotype=phenotype, "locus ID"=private$data$locus$name()))},

		## match partial names, check and throw error for ambiguous and unknown, remove duplicates/expand wildcards
		check.phenotypes = function(phenotypes) {
			if (length(phenotypes) > 0) return(check.phenotypes(phenotypes, names(private$phenotype.info$status), label="processed locus"))
			else if (is.null(phenotypes)) return(NULL)
			else return(character(0))
		},

# NB: currently using trim.components on component drop, which changes shared ld.block object (TODO)
		process = function(locus.data) {
			phenotypes = locus.data$phenotypes(); estimates = list()
			status = named.list(phenotypes, rep("available", length(phenotypes)))

			while (TRUE) { ## rerun estimators for all (non-failed) phenotypes until all completed or failed
				rerun = FALSE
				for (ph in phenotypes[unlist(status) != "failed"]) {
					curr.stats = locus.data$get.sumstats()[[ph]]
					res = tryCatch(
						curr.stats$compute.marginal(),
						dropped.components = function(err) {return(err)},
						data.error = function(err) {private$warning(err, phenotype=ph); return(err)},
						error = function(err) {private$failure(err, phenotype=ph)}
					)
					if (is.error(res)) {
						if (is.error(res, "dropped.components")) {   ## break out of for-loop, goes to next iteration of while loop and reruns all (non-failed) phenotypes
							ld.block = locus.data$get.ld()$trim.components(drop=res$get.data()$dropped) ## NB: currently, this modified the ld.block IN PLACE
							if (ld.block$no.components() < private$settings$get("minimum.components")) private$failure("fewer than ", private$settings$get("minimum.components"), " components in locus after pruning unstable components")
							rerun = TRUE; break
						}
						if (is.error(res, "failure")) status[[ph]] = "failed"
					}	else estimates[[ph]] = res
				}
				if (!rerun) break  ## all phenotypes processed or failed, break out of while loop
			}

			estimates = estimates[phenotypes[unlist(status) != "failed"]]
			neg.sigma = as.logical(unlist(lapply(estimates, function(est) {is.na(est$get("sigma")) || est$get("sigma") <= 0})))
			if (any(neg.sigma)) private$warning("%phenotype% ", items.and=names(estimates)[neg.sigma], " %has% a negative or zero sampling variance")
			status[names(status) %in% names(estimates)[neg.sigma]] = "failed"

			neg.omega = !neg.sigma & unlist(lapply(estimates, function(est) {is.na(est$get("omega")) || est$get("omega") <= 0}))
			if (any(neg.omega)) private$warning("%phenotype% ", items.and=names(estimates)[neg.omega], " %has% a negative or zero genetic variance estimate")
			status[names(status) %in% names(estimates)[neg.omega]] = "negative"

			estimates = estimates[phenotypes[unlist(status) != "failed"]]
			if (length(phenotypes) > 0 && length(estimates) == 0) private$failure("processing has failed for all phenotypes")

			private$estimates = estimates
			private$phenotype.info$status = status
			private$data$correlations = locus.data$get.correlations()$subset(names(estimates))
		},

		subset.internal = function(phenotypes) {
			for (param in names(private$phenotype.info)) private$phenotype.info[[param]] = private$phenotype.info[[param]][phenotypes]

			pheno.processed = phenotypes[phenotypes %in% names(private$estimates)]
			private$estimates = private$estimates[pheno.processed]
			private$data$correlations = private$data$correlations$subset(pheno.processed)
			return(self)
		},

		composite.internal = function(names, phenotypes, weights, discard.input=F) {
			for (i in seq_along(names)) {
				curr.name = names[i]
				input = private$estimates[phenotypes[[i]]]

				if (weights[[i]][1] == "meta")	estimate = MetaAnalysisEstimates$new(names[i], input, private$data$correlations)
				else estimate = CompositeEstimates$new(names[i], input, private$data$correlations, weights[[i]])

				private$phenotype.info$info[[curr.name]] = estimate$phenotype.info()
				if (!is.na(estimate$get("sigma")) && estimate$get("sigma") > 0) {
					private$estimates[[curr.name]] = estimate
					if (!is.na(estimate$get("omega")) && estimate$get("omega") > 0) {
						private$phenotype.info$status[[curr.name]] = "available"
						private$data$correlations = estimate$get.correlations()
					} else {
						private$phenotype.info$status[[curr.name]] = "negative"
						private$warning("composite phenotype '", curr.name, "' has a negative or zero genetic variance estimate")
					}
				} else {
					private$phenotype.info$status[[curr.name]] = "failed"
					private$warning("composite phenotype '", curr.name, "' has a negative or zero sampling variance")
				}
			}

			if (discard.input) {
				keep = names(private$phenotype.info$status)[!(names(private$phenotype.info$status) %in% unique(unlist(phenotypes)))]
				private$subset.internal(keep)
			}
			return(self)
		}
	),
	public = list(
		initialize = function(locus.data, ...) {
			private$settings = initialize.settings(self, ...)
			check.types(locus.data="LocusData")

			private$phenotype.info$info = lapply(locus.data$get.sumstats(), function(ss) {ss$phenotype.info()})
			private$data$locus = locus.data$get.locus()
			private$data$ld.block = locus.data$get.ld()

			private$process(locus.data)
		},

		print = function(...) {
			summary = self$summarize()
			printer = ObjectPrinter$new("ProcessedLocus")$parameter.list(snps=private$data$ld.block$no.snps(), components=private$data$ld.block$no.components())
			printer$add.parameter("phenotypes", length(private$phenotype.info$status))
			if (!is.null(summary) && nrow(summary) > 0) {
				for (i in 1:nrow(summary)) {
					phenotype = summary$phenotype[i]
					type = summary$type[i]
					status = if (summary$status[i] != "available") toupper(summary$status[i])
					printer$add.line(phenotype, parentheses(type), parentheses(status), indent=2)
				}
			}
			cat(printer$to.string())
		},

		phenotypes = function(include=c("all", "processed", "available", "positive", "failed")) {
			include = match.arg(include); phenotypes = as.character(names(private$phenotype.info$status))
			if (length(phenotypes) > 0) {
				if (include == "all") return(phenotypes)
				if (include == "processed") return(phenotypes[unlist(private$phenotype.info$status) != "failed"])
				if (include == "positive" || include == "available") return(phenotypes[unlist(private$phenotype.info$status) == "available"])
				if (include == "failed") return(phenotypes[unlist(private$phenotype.info$status) == "failed"])
			} else return(character(0))
		},

		## return list of matched phenotype names (check validity, convert partial names to full, remove duplicates, expand wildcards)
		match.phenotypes = function(phenotypes) {return(if (length(phenotypes) > 0) private$check.phenotypes(phenotypes) else character(0))},

		no.pheno = function(include=c("all", "processed", "available", "positive", "failed")) {return(length(self$phenotypes(include)))},
		no.snps = function() {return(private$data$ld.block$no.snps())},
		no.components = function() {return(private$data$ld.block$no.components())},

		get.locus = function() {return(private$data$locus)},
		get.ld = function() {return(private$data$ld.block)},

		get.correlations = function(incl.negative=F) {
			if (!incl.negative) return(private$data$correlations$subset(self$phenotypes("positive")))
			else return(private$data$correlations)
		},
		get.estimates = function(incl.negative=F) {return(if (!incl.negative) private$estimates[self$phenotypes("positive")] else private$estimates)},

		## provide overview of status and marginal statistics
		summarize = function() {
			if (length(private$phenotype.info$status) > 0) {
				summary = data.frame(
					locus = private$data$locus$name(),
					phenotype = names(private$phenotype.info$status),
					type = unlist(lapply(private$phenotype.info$info, function(pi) {pi$get.traits(sep=";")})),
					status = unlist(private$phenotype.info$status),
					no.snps = NA, no.components = NA,
					h2.observed=NA, h2.latent=NA, h2.composite=NA, p.univariate=NA,
					stringsAsFactors=F
				)

				index = match(names(private$estimates), summary$phenotype)
				summary$h2.observed[index] = unlist(lapply(private$estimates, function(est) {est$get("h2")}))
				summary$p.univariate[index] = unlist(lapply(private$estimates, function(est) {est$get("p.value")}))
				summary$no.snps[index] = unlist(lapply(private$estimates, function(est) {est$no.snps()}))
				summary$no.components[index] = unlist(lapply(private$estimates, function(est) {est$no.components()}))

				if (any(sapply(private$phenotype.info$info, function(pi) {pi$is("binary")}))) summary$h2.latent[index] = unlist(lapply(private$estimates, function(est) {est$get("h2.latent")}))
				else summary$h2.latent = NULL

				if (any(sapply(private$phenotype.info$info, function(pi) {pi$is("composite") && !pi$is("meta")}))) summary$h2.composite[index] = unlist(lapply(private$estimates, function(est) {est$get("h2.composite")}))
				else summary$h2.composite = NULL

				return(add.rownames(summary, NULL))
			} else return(NULL)
		},

		## keep only specified phenotypes if include=T, discard specified phenotypes otherwise
		subset = function(phenotypes, include=T, empty.mode=c("warning", "error", "quiet")) {empty.mode = match.arg(empty.mode)
			if (length(phenotypes) == 0) phenotypes = character(0)
			phenotypes = private$check.phenotypes(phenotypes)
			if (!include) phenotypes = names(private$phenotype.info$status)[!(names(private$phenotype.info$status) %in% phenotypes)]
			if (length(phenotypes) == 0 && self$no.pheno() > 0 && empty.mode != "quiet") ifelse(empty.mode == "warning", private$warning, private$failure)("no phenotypes remaining after subsetting")

			return(self$clone()$.__enclos_env__$private$subset.internal(phenotypes))
		},

		filter = function(include=c("available", "positive", "processed", "failed"), empty.mode=c("warning", "error", "quiet")) {return(self$subset(self$phenotypes(match.arg(include)), include=T, empty.mode=match.arg(empty.mode)))},
		filter.h2 = function(min.h2, empty.mode=c("warning", "error", "quiet")) {
			summary = self$summarize(); phenotypes = c()
			if (!is.null(summary)) {
				top.h2 = apply(summary[grep("^h2", names(summary))], 1, function(h2) {ifelse(all(is.na(h2)), NA, max(h2, na.rm=T))})
				phenotypes = summary$phenotype[top.h2 & top.h2 >= min.h2]
			}
			return(self$subset(phenotypes, include=T, empty.mode=match.arg(empty.mode)))
		},
		filter.pval = function(max.pval, empty.mode=c("warning", "error", "quiet")) {
			summary = self$summarize(); phenotypes = c()
			if (!is.null(summary)) phenotypes = summary$phenotype[!is.na(summary$p.univariate) & summary$p.univariate <= max.pval]
			return(self$subset(phenotypes, include=T, empty.mode=match.arg(empty.mode)))
		},

		## 'phenotypes' should be a vector of phenotypes (does partial matching)
		## 'weights' either be NULL/1, "meta"/"ivw" (both are equivalent), or a numeric vector with the same length as 'phenotypes'
		##  - negative weights are allowed, but this disables heritability estimates
		## to construct multiple composite phenotypes at once, supply a vector of names, and list arguments for 'phenotype' and 'weights' of equal length
		## if discard.input=T, phenotypes used to contruct the composites are removed from the ProcessedLocus object
		add.composite = function(names, phenotypes, weights=NULL, discard.input=F) {
			if (length(names) == 0) private$error("no %name% provided for composite %phenotype%", .plural=is.list(phenotypes) && length(phenotypes) > 1)
			if (any(duplicated=duplicated(names))) private$error("'names' argument contains duplicate values")
			IF.ANY(exist=names %in% names(private$phenotype.info$status), THEN=private$error("%name% ", items.and=names[exist], " %is% already in use"))

			if (!is.list(phenotypes)) phenotypes = list(phenotypes)
			if (length(phenotypes) != length(names)) private$error("number of phenotype vectors does not match number of composite phenotype names")
			phenotypes = lapply(phenotypes, private$check.phenotypes)
			IF.ANY(unavailable = !(unlist(phenotypes) %in% self$phenotypes("available")), THEN=private$error("%phenotype% ", items.and=unique(unlist(phenotypes)[unavailable]), " used for composite phenotypes %is% not available"))

			if (!is.null(weights)) {
				if (!is.list(weights)) weights = if (length(names) == 1) list(weights) else as.list(weights)
				if (length(weights) != length(names)) private$error("number of weight specifications does not match number of requested composite phenotypes")
				weights = lapply(weights, function(wt) {if (any(c("meta", "ivw") %in% tolower(wt))) "meta" else wt})
				no.weights = sapply(weights, length)
				IF.ANY(mismatch = !(no.weights == sapply(phenotypes, length) | no.weights == 1), THEN=private$error("number of weights provided for composite %phenotype% ", items.and=names[mismatch], " does not match corresponding number of input phenotypes"))
			} else weights = as.list(rep(1, length(names)))

			return(self$clone()$.__enclos_env__$private$composite.internal(names, phenotypes, weights, discard.input))
		}
	)
)











