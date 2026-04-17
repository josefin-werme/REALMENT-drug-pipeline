#' @include global.R
#' @include model__core.R
#' @include model__bivariate.R


MultivariateAnalysis = R6::R6Class("MultivariateAnalysis",
	inherit = AnalysisModel,
	private = list(
		settings = DEFAULT.SETTINGS("MultivariateAnalysis",
			coefficient.limit = add.info(globals$regression.invalid.bounds, type="numeric", min=1),   ## standardized regression coefficient bound (considered invalid if exceeded)
			maximum.r2 = add.info(globals$partial.cor.max.r2, type="numeric", range=c(0,0.99))   ## maximum explained variance for univariate models in partial correlation analysis
		),

		get.input = function(outcomes, predictors) {
			input = list(no.components = private$data$no.components(), omega = private$data$get.omega(outcomes, predictors), sigma = private$data$get.sigma(outcomes, predictors))
			input$omega.x = input$omega[predictors,predictors,drop=F]; input$omega.xy = input$omega[predictors,outcomes,drop=F]; input$omega.y = input$omega[outcomes,outcomes,drop=T]
			input$omega.x.inv = muffleErrors(eigen.invert(input$omega.x))
			input$rg = private$data$get.rg(outcomes, predictors, missing.bound=private$settings$get("rg.limit"))

			info = list(
				model.name = paste0(paste(outcomes, collapse=" ~ "), " | ", paste(predictors, collapse=" + ")),
				outcome = outcomes, predictors = predictors,
				no.snps = private$data$no.snps(), no.components = input$no.components
			)

			status = list()
			pair.types = unique(as.vector(private$data$get.types(outcomes, predictors)))
			unavailable = grep("^unavailable", pair.types, value=T)
			if (length(unavailable) > 0) status$model = paste(gsub("unavailable - ", "", unavailable), collapse="; ")
			else if ("truncated" %in% pair.types) status$model = "model contains truncated phenotypes"
			else if (any(is.na(input$omega))) status$estimation = "missing values in genetic covariance matrix"
			else if (any(is.na(input$rg[lower.tri(input$rg)]))) status$estimation = "genetic correlation matrix contain out of bounds values"
			else if (is.error(input$omega.x.inv)) status$estimation = "predictor genetic covariance matrix is not invertible"

			return(list(input=input, info=info, status=status))
		}
	)
)



## outcome should be a single variable; outcome is removed from predictor set if present
MultipleRegression = R6::R6Class("MultipleRegression",
	inherit = MultivariateAnalysis,
	private = list(
		compute.model = function(outcome, predictors) {
			data = private$get.input(outcome, predictors); input = data$input; status = data$status

			coefficients = data.frame(outcome=outcome, predictor=predictors, gamma.raw=NA, gamma=NA)
			parameters = data.frame(data$info[c("model.name", "outcome")], predictors=I(list(predictors)), no.predictors=length(predictors), data$info[c("no.components", "no.snps")])
			parameters[c("tau.raw", "tau", "r2")] = NA
			if (private$settings$get("compute.ci") != "never") {coefficients[c("gamma.lower", "gamma.upper")] = NA; parameters[c("r2.lower", "r2.upper")] = NA}
			if (private$settings$get("compute.pval")) {coefficients$test.type = "sampling"; coefficients[c("sampling.variance", "no.iter", "p.value")] = NA}

			if (length(status) == 0) {
				if (length(predictors) > 1) {
					coefficients$gamma.raw = as.vector(input$omega.x.inv %*% input$omega.xy)
					coefficients$gamma = coefficients$gamma.raw * sqrt(diag(input$omega.x)/input$omega.y)

					coefficients$gamma[abs(coefficients$gamma) > private$settings$get("coefficient.limit")] = NA
					if (!any(is.na(coefficients$gamma))) {
						parameters$tau.raw = as.numeric(input$omega.y - t(input$omega.xy) %*% input$omega.x.inv %*% input$omega.xy)
						parameters$tau = parameters$tau.raw / input$omega.y
						parameters$r2 = 1 - parameters$tau

						phenotypes = c(predictors, outcome)
						if (private$settings$get("compute.pval")) {
							res = muffleErrors(integral.p(multivariate.integral, K=input$no.components, omega=input$omega[phenotypes,phenotypes]/input$no.components, sigma=input$sigma[phenotypes,phenotypes], adap.thresh=private$settings$get("adaptive.threshold")))
							if (!is.error(res)) coefficients[c("no.iter", "p.value")] = res[c("no.iter", "p.value")]
							else status$p.value = error.message(res)
						}

						if (private$settings$get("compute.ci") == "significant") compute.ci = !all(is.na(coefficients$p.value)) & any(coefficients$p.value <= private$settings$get("ci.threshold"))
						else compute.ci = (private$settings$get("compute.ci") == "always")
						if (compute.ci) {
							res = muffleErrors(ci.multivariate(K=input$no.components, omega=input$omega[phenotypes,phenotypes]/input$no.components, sigma=input$sigma[phenotypes,phenotypes]))
							if (!is.error(res)) {
								coefficients$gamma.lower = res$gamma.lower
								coefficients$gamma.upper = res$gamma.upper
								parameters[c("r2.lower", "r2.upper")] = res$r2.interval
							} else status$ci = error.message(res)
						}
					} else status$estimation = "coefficient estimates are out of bounds"
				} else {
					bivariate = BivariateCorrelation$new(private$data, private$settings, rg.limit=private$settings$get("coefficient.limit"))$compute(outcome, predictors)
					status = bivariate$status

					coefficients$gamma.raw = input$omega[predictors,outcome] / input$omega[predictors,predictors]
					coefficients$gamma = bivariate$results$rho
					parameters$tau.raw = input$omega[outcome,outcome] - input$omega[predictors,outcome]^2 / input$omega[predictors,predictors]
					parameters$tau = 1 - bivariate$results$r2
					parameters$r2 = bivariate$results$r2

					if (private$settings$get("compute.pval")) {
						cols = c("test.type", "sampling.variance", "no.iter", "p.value"); cols = cols[cols %in% names(bivariate$results)]
						coefficients[cols] = bivariate$results[cols]
					}

					if (private$settings$get("compute.ci") != "never") {
						coefficients[c("gamma.lower", "gamma.upper")] = bivariate$results[c("rho.lower", "rho.upper")]
						parameters[c("r2.lower", "r2.upper")] = bivariate$results[c("r2.lower", "r2.upper")]
					}
				}
			}

			if (all(coefficients$test.type == "sampling")) coefficients$sampling.variance = NULL
			else if (!any(coefficients$test.type == "sampling")) coefficients$no.iter = NULL

			return(list(results=list(parameters=parameters, coefficients=coefficients), status=status))
		}
	),
	public = list(
		compute = function(outcome, predictors) {
			outcome = private$resolve.phenotypes(outcome)
			if (length(outcome) != 1) input.error(ifelse(length(outcome) > 1, "multiple", "no"), " outcome phenotypes specified for regression model")
			predictors = sort(private$resolve.phenotypes(predictors, exclude=outcome))
			if (length(predictors) == 0) input.error("no predictor phenotypes specified for regression model")

			return(private$compute.model(outcome, predictors))
		}
	)
)


## outcomes can contain more than two phenotypes, will run all pairs of outcomes if so
## outcomes are filtered out of the predictor set for analyses they are involved in only, not for other those of other outcome pairs
PartialCorrelation = R6::R6Class("PartialCorrelation",
	inherit = MultivariateAnalysis,
	private = list(
		marginal = list(),

		compute.model = function(outcomes, predictors) {
			data = private$get.input(outcomes, predictors); input = data$input; status = data$status

			output = data.frame(
				outcome1 = outcomes[1],	outcome2 = outcomes[2],
				predictors = I(list(predictors)),
				pair.type = private$data$get.types(outcomes)[1,2],
				no.predictors = length(predictors),
				data$info[c("no.snps", "no.components")],
				r2.marginal1=NA, r2.marginal2=NA,
				rho=data$input$rg[outcomes[1],outcomes[2]], rho.partial=NA,
				stringsAsFactors=F
			)

			if (private$settings$get("compute.ci") != "never") output[c("rho.partial.lower", "rho.partial.upper")] = NA
			if (private$settings$get("compute.pval")) output[c("no.iter", "p.value")] = NA

			if (length(status) == 0) {
				if (abs(output$rho) > private$settings$get("rg.limit")) output$rho = NA
				if (is.na(output$rho)) status$estimation = "marginal correlation is out of bounds"
			}

			if (length(status) == 0) {
				labels = paste0(outcomes, " | ", paste(predictors, collapse=" + "))
				for (i in 1:2) {
					if (!(labels[i] %in% names(private$marginal))) {
						pheno = outcomes[i]
						marginal = MultipleRegression$new(private$data, private$settings, coefficient.limit=Inf, compute.pval=F, compute.ci="never")$compute(pheno, predictors)
						if (is.null(marginal$status$estimation)) {
							if (marginal$results$parameters$r2 <= private$settings$get("maximum.r2")) {
								marginal$partial.variance = input$omega.y[pheno,pheno] - t(input$omega.xy[,pheno]) %*% input$omega.x.inv %*% input$omega.xy[,pheno]
								if (is.na(marginal$partial.variance) || marginal$partial.variance <= 0) marginal$status$estimation = "invalid partial variance"
							} else marginal$status$estimation = paste0("marginal r2 exceeds maximum of ", private$settings$get("maximum.r2"))
						} else marginal$status$estimation = "unable to compute marginal r2"
						private$marginal[[labels[i]]] = marginal
					}
				}


				univ.r2 = sapply(private$marginal[labels], function(m) {ifelse(length(m$status) == 0, m$results$parameters$r2, NA)})
				if (!any(is.na(univ.r2))) {
					output[,c("r2.marginal1", "r2.marginal2")] = univ.r2
					partial.variance = list.extract(private$marginal[labels], "partial.variance", flatten=T)
					partial.covariance = input$omega.y[outcomes[1],outcomes[2]] - t(input$omega.xy[,outcomes[1]]) %*% input$omega.x.inv %*% input$omega.xy[,outcomes[2]]
					output$rho.partial = as.numeric(partial.covariance / sqrt(prod(partial.variance)))
					if (abs(output$rho.partial) > private$settings$get("rg.limit")) output$rho.partial = NA

					if (!is.na(output$rho.partial)) {
						phenotypes = c(predictors, outcomes)
						if (private$settings$get("compute.pval")) {
							res = muffleErrors(integral.p(pcov.integral, K=input$no.components, omega=input$omega[phenotypes,phenotypes]/input$no.components, sigma=input$sigma[phenotypes,phenotypes], adap.thresh=private$settings$get("adaptive.threshold")))
							if (!is.error(res)) output[c("no.iter", "p.value")] = res[c("no.iter", "p.value")]
							else status$p.value = error.message(res)
						}

						if (private$settings$get("compute.ci") == "always" || (private$settings$get("compute.ci") == "significant" && !is.na(output$p.value) && output$p.value <= private$settings$get("ci.threshold"))) {
							res = muffleErrors(ci.pcor(K=input$no.components, omega=input$omega[phenotypes,phenotypes]/input$no.components, sigma=input$sigma[phenotypes,phenotypes]))
							if (!is.error(res)) output[c("rho.partial.lower", "rho.partial.upper")] = res
							else status$ci = error.message(res)
						}
					} else status$estimation = "partial correlation is out of bounds"
				} else status$estimation = "invalid marginal model"
			}

			return(list(info = data$info, output=output, status=status))
		}
	),
	public = list(
		compute = function(outcomes, predictors) {
			outcomes = sort(private$resolve.phenotypes(outcomes))
			if (length(outcomes) < 2) input.error(ifelse(length(outcomes) == 1, "one", "no"), " outcome %phenotype% specified for partial correlation model", .plural=length(outcomes) != 1)
			predictors = sort(private$resolve.phenotypes(predictors, exclude=if (length(outcomes) == 2) outcomes))
			if (length(predictors) == 0) input.error("no predictor phenotypes specified for partial correlation model")

			results = list(); status = list(); private$marginal = list()
			for (phs in make.pairs(outcomes, allow.equal=F, ignore.order=T, sort=T)) {
				if (!all(predictors %in% phs)) {
					curr = private$compute.model(phs, predictors[!(predictors %in% phs)])
					label = curr$info$model.name
					results[[label]] = data.frame(model.name=label, curr$output, stringsAsFactors=F)
					if (length(curr$status) > 0) status[[label]] = curr$status
				}
			}
			results = add.rownames(fill.rowmerge(results), names=NULL)

			marginal.status = list.extract(private$marginal, "status", trim=T)
			return(list(results=results, status=status, marginal.status=marginal.status))
		}
	)
)






