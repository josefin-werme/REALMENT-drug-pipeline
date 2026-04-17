#' @include global.R
#' @include processing__locus.R
#' @include model__core.R
#' @include model__bivariate.R
#' @include model__multivariate.R


## general analysis interface object, exposes all primary analysis models
## NB: setting ci.threshold explicitly will automatically switch compute.ci="significant"

## STATUS: volatile, deep (ProcessedLocus, LocusModel)
AnalysisProcessor = R6::R6Class("AnalysisProcessor",
	private = list(
		settings = DEFAULT.SETTINGS("AnalysisProcessor",
			IMPORT.DEFAULTS("AnalysisModel", "MultivariateAnalysis"),
			univ.threshold = add.info(1, type="numeric", range=c(0,1)),  ## univariate p-value threshold for phenotype inclusion
			signif.decimals = globals$analysis.signif.decimals,   ## number of significant decimals to round estimates off to
			extra.columns = c("info", "none", "statistical", "all"),  ## which additional columns to include, as listed below; "statistical" includes info.columns as well
			info.columns = add.info(c("no.phenotypes", "no.snps", "no.components"), type="character", vector=T),
			statistical.columns = add.info(c("omega.*", "tau.*", "*.raw", "sampling.variance", "test.type", "no.iter"), type="character", vector=T),
			technical.columns = add.info(c("id", "pair.type"), type="character", vector=T)
		),

		data = NULL,   ## ProcessedLocus object
		model = NULL,   ## LocusModel object
		univ.info = NULL,   ## data.frame obtained data$summarize(); status is updated to selected for available phenotypes with p-value below univ.threshold

		## if lines is a named list, will print "[name] ([value])"
		format.condition = function(type.func, ..., phenotype=NULL, lines=NULL) {
			if (length(lines) > 0 && is.list(lines)) {
				values = unlist(lines)
				lines = paste0(names(lines), ifelse(nchar(values) > 0, paste0(" (", values, ")"), values))
			}
			label = list(phenotype=phenotype, "locus ID"=private$data$get.locus()$name())
			type.func(..., .label=label, .lines=lines)
		},

		warning = function(..., phenotype=NULL, lines.list=NULL) {private$format.condition(throw.warning, ..., phenotype=phenotype, lines=lines.list)},
		failure = function(..., phenotype=NULL, lines.list=NULL) {private$format.condition(data.error, ..., phenotype=phenotype, lines=lines.list)},
		error = function(..., phenotype=NULL, lines.list=NULL) {private$format.condition(input.error, ..., phenotype=phenotype, lines=lines.list)},


		## check if at least two phenotypes in data
		check.size = function(abort=F) {
			if (private$model$no.pheno() <= 1) {
				msg = list("%[?only one/no]% %phenotype% available in input data", if (private$settings$get("univ.threshold") < 1) " at univariate p-value threshold of ", private$settings$get("univ.threshold"), if (abort) ", aborting analysis", .plural=(private$model$no.pheno() == 0))
				ifelse(abort, private$failure, private$warning)(msg)
			}
		},

		## check validity of phenotypes, match partial phenotype names against available phenotypes, and remove duplicates/expand wildcards; NB: NULL input returns character(0) output
		check.phenotypes = function(phenotypes) {
			if (length(phenotypes) > 0) return(check.phenotypes(phenotypes, private$univ.info$phenotyp, label="analysis input"))
			else return(character(0))
		},

		## check validity of phenotypes, then check which are available for analysis in locus.model
		## if error.none=T, will throw an error if all are missing
		## if missing="quiet", missing phenotypes will be dropped quietly, otherwise a warning/error is also generated
		check.available = function(phenotypes, label=NULL, missing=c("quiet", "warning", "error"), error.none=T) {missing = match.arg(missing)
			if (length(phenotypes) > 0) {
				phenotypes = private$check.phenotypes(phenotypes)
				available = phenotypes %in% private$model$phenotypes()
				if (!any(available) && (error.none || missing == "error")) private$failure("%[+none of]% specified ", label, " %phenotype% %is% %[-not]% available for analysis", .plural=sum(!available))
				if (!all(available) && missing != "quiet")	ifelse(missing == "warning", private$warning, private$failure)(label, " %phenotype% ", items.and=phenotypes[!available], " %is% not available for analysis")
				return(phenotypes[available])
			} else return(character(0))
		},

		set.input = function() {
			if (private$data$no.pheno() > 0) {
				private$univ.info = private$data$summarize()
				private$univ.info$status[private$univ.info$status == "available" & (is.na(private$univ.info$p.univariate) | private$univ.info$p.univariate > private$settings$get("univ.threshold"))] = "non-significant"
			} else private$univ.info = data.frame(phenotype=character(0), status=character(0))

			private$model = LocusModel$new(private$data$subset(self$phenotypes("available"), empty.mode="quiet"))
			private$check.size()
		},

		process.settings = function(...) {
			settings = if (is.null(private$settings)) initialize.settings(self, ...) else private$settings$set(..., .in.place=F)
			args = list(...)
			if (!is.null(names(args)) && "ci.threshold" %in% complete.names(names(args), names(settings$get.all()))) settings$set(compute.ci = "significant")
			return(settings)
		},

		## expects named (model.name) list, with named list of failure (type=reason)
		## types "model" and "estimation" are considered failures, other types are considered warnings
		## returns vector of names for failed models
		parse.status = function(status.list, model.type, is.plural) {
			if (length(status.list) > 0) {
				failures = trim.list(lapply(status.list, function(st) {if (any(names(st) %in% c("model", "estimation"))) st[names(st) %in% c("model", "estimation")][1]}))
				if (length(failures) > 0) {
					failures = lapply(failures, function(f) {paste0(names(f)[1], " failed: ", f[[1]])})
					private$warning("unable to compute following ", model.type, " %model%", .plural=length(failures), lines.list=failures)
				}

				incomplete = trim.list(status.list[!(names(status.list) %in% names(failures))])
				if (length(incomplete) > 0) {
					for (type in c("ci", "p.value")) {
						has.type = unlist(lapply(incomplete, function(s) {length(s[[type]]) > 0}))
						if (any(has.type)) {
							label = switch(type, estimation="estimate", ci="confidence interval", p.value="p-value")
							private$warning(pluralize("unable to compute %[", label, "]%", .plural=is.plural), " for following ", model.type, " %model%", .plural=sum(has.type), lines.list=incomplete[has.type])
						}
					}
				}
				return(names(failures))
			}
		},

		## filter columns, and set named arguments with 'column=bounds', where bounds is considered symmetric around zero if only a single value (ignores columns if not present)
		## truncates columns with given name to specified bounds and rounds off to settings$signif.decimals with signif() (set bound to Inf to round only)
		## applies bounds to columns [NAME].[SUFFIX] as well, unless column name is explicitly set to NA (eg. rho=1 will apply bounds to rho.lower as well, unless rho.lower=NA argument is present)
		format.results = function(results, ..., settings=NULL) {
			if (is.null(settings)) settings = private$settings
			params = lapply(list(...), function(x) {return(if(length(x) == 1) c(-x,x) else x[1:2])})
			columns = params[names(params) %in% names(results)]
			for (par in names(params)[order(nchar(names(params)), decreasing=T)]) { ## do in descending order of length, so best match wins
				matches = grep(paste0("^", par, "\\."), names(results), value=T)
				for (m in matches[!(matches %in% names(columns))]) columns[[m]] = params[[par]]
			}

			for (col in names(columns[!unlist(lapply(columns, function(c) {any(is.na(c))}))])) {
				curr = results[[col]]
				curr[!is.na(curr) & curr < columns[[col]][1]] = columns[[col]][1]
				curr[!is.na(curr) & curr > columns[[col]][2]] = columns[[col]][2]
				results[[col]] = signif(curr, settings$get("signif.decimals"))
			}

			filter = c()
			if (settings$get("extra.columns") != "all") filter = c(filter, settings$get("technical.columns"))
			if (!(settings$get("extra.columns") %in% c("all", "statistical"))) filter = c(filter, settings$get("statistical.columns"))
			if (settings$get("extra.columns") == "none") filter = c(filter, settings$get("info.columns"))

			if (length(filter) > 0) {
				drop = names(results) %in% filter
				for (wildcard in grep("*", filter, value=T, fixed=T)) drop = drop | grepl(wildcard, names(results))
				results = results[!drop]
			}

			return(add.rownames(results, names=NULL))
		}
	),
	public = list(
		initialize = function(data, ...) {
			check.types(data="ProcessedLocus")
			private$settings = private$process.settings(...)

			private$data = data
			private$set.input()
		},

		no.pheno = function(include=c("all", "available")) {return(length(self$phenotypes(match.arg(include))))},
		phenotypes = function(include=c("all", "available")) {include = match.arg(include)
			if (include == "available") return(private$univ.info$phenotype[private$univ.info$status == "available"])
			else return(private$univ.info$phenotype)
		},

		get.settings = function(as.list=F) {return(if (as.list) private$settings$get.all() else private$settings$clone())},

		get.locus.data = function() {return(private$data)},
		get.locus.model = function() {return(private$model)},

		get.univariate = function(phenotypes=NULL) {
			if (length(phenotypes) > 0) return(add.rownames(private$univ.info[match(private$check.phenotypes(phenotypes), private$univ.info$phenotype),], names=NULL))
			else return(private$univ.info)
		},

		get.rg = function(phenotypes=NULL, truncate=T, rg.limit=NULL) {
			if (!is.null(phenotypes)) phenotypes = private$check.available(phenotypes, missing="warning")
			return(private$model$get.rg(phenotypes, truncate=truncate, missing.bound=private$settings$validate("rg.limit", rg.limit)))
		},

		set.parameter = function(...) {
			univ.threshold = settings$get("univ.threshold")
			private$settings = private$process.settings(...)
			if (private$settings$get("univ.threshold") != univ.threshold) private$set.input()
			invisible(self)
		},


		### Main Analysis Functions ###
		##  NB: settings arguments in ... are only used for requested analysis and cannot change the AnalysisProcessor object itself
		##  ie. 'univ.threshold' argument will be ignored, to change threshold call set.parameter() first instead

		## runs correlation for all phenotypes in 'targets' with those in 'phenotypes'
		## if targets=NULL it defaults to be equal to 'phenotypes' (ie. all phenotypes if phenotypes=NULL)
		## if phenotypes=NULL it defaults to all available phenotypes
		bivariate = function(targets=NULL, phenotypes=NULL, ...) {settings = private$process.settings(...)
			private$check.size(abort=T); private$check.available(c(phenotypes, targets), missing="warning")

			phenotypes = if (!is.null(phenotypes)) private$check.available(phenotypes) else self$phenotypes("available")
			targets = if (!is.null(targets)) private$check.available(targets, label="target") else phenotypes
			if (length(phenotypes) == 1 && length(targets) == 1 && phenotypes == targets) private$error("no phenotype pairs to analyze")

			output = BivariateCorrelation$new(private$model, settings)$compute(targets=targets, phenotypes=phenotypes)
			failed = private$parse.status(output$status, model.type="bivariate correlation", is.plural=F)
			return(private$format.results(output$results, rho=1, r2=c(0,1), standard.error=Inf, p.value=Inf, settings=settings)[!(output$results$model.name %in% failed),])
		},

		## runs multiple regressions of each phenotype in 'outcomes' (separately) on 'predictors'
		## if predictors=NULL it will default to all available phenotypes
		## for any outcome, if it is contained in predictors it will be removed from predictors for that analysis, but not when analysing other outcomes
		regression = function(outcomes, predictors=NULL, ...) {settings = private$process.settings(...)
			private$check.size(abort=T)

			outcomes = private$check.available(outcomes, label="outcome", missing="warning")
			predictors = if (!is.null(predictors)) private$check.available(predictors, label="predictor", missing="warning") else self$phenotypes("available")

			if (length(predictors) == 1 && any(outcomes %in% predictors)) {
				if (length(outcomes) == 1) private$failure("no valid predictor phenotypes have been specified")
				else private$warning("no valid predictor phenotypes available for outcome '", predictors, "'")
				outcomes = outcomes[outcomes != predictors]
			}

			output = list.extract(lapply(outcomes, function(o) {MultipleRegression$new(private$model, settings)$compute(outcome=o, predictors=predictors[predictors != o])}), "results$parameters", "results$coefficients", "status")
			model.names = list.extract(output$parameters, "model.name", flatten=T)
			failed = private$parse.status(add.names(output$status, names=model.names), model.type="regression", is.plural=T)

			parameters = fill.rowmerge(output$parameters)
			coefficients = fill.rowmerge(lapply(seq_along(model.names), function(i) {cbind(model.name=model.names[i], output$coefficients[[i]])}))

			results = list(
				models = private$format.results(parameters[!(parameters$model.name %in% failed),], tau.raw=c(0,Inf), tau=c(0,1), r2=c(0,1), settings=settings),
				coefficients = private$format.results(coefficients[!(coefficients$model.name %in% failed),], gamma.raw=Inf, gamma=settings$get("coefficient.limit"), p.value=Inf, settings=settings)
			)

			return(results)
		},


		## runs partial correlation for all pairs of phenotypes in 'outcomes'
		## if predictors=NULL it will default to all available phenotypes
		## for any outcome, if it is contained in predictors it will be removed from predictors for every analysis it is part of, but not when analysing other outcomes
		partial = function(outcomes, predictors=NULL, ...) {settings = private$process.settings(...)
			private$check.size(abort=T)

			outcomes = private$check.available(outcomes, label="outcome", missing="warning")
			if (length(outcomes) == 1) private$failure("only one valid outcome phenotype specified")

			predictors = if (!is.null(predictors)) private$check.available(predictors, label="predictor", missing="warning") else self$phenotypes("available")
			if (length(outcomes) == 2) predictors = predictors[!(predictors %in% outcomes)]
			if (length(predictors) == 0) private$failure("no valid predictor phenotypes have been specified")

			output = PartialCorrelation$new(private$model, settings)$compute(outcomes=outcomes, predictors=predictors)
			failed = private$parse.status(output$status, model.type="partial correlation", is.plural=F)
			return(private$format.results(output$results, rho=1, r2=c(0,1), p.value=Inf, settings=settings)[!(output$results$model.name %in% failed),])
		}
	)
)



