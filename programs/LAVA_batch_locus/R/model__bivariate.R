#' @include global.R
#' @include model__core.R



BivariateCorrelation = R6::R6Class("BivariateCorrelation",
	inherit = AnalysisModel,
	private = list(
		compute.rg = function(pheno1, pheno2) {
			no.components = private$data$no.components()
			omega = private$data$get.omega(pheno1, pheno2)
			sigma = private$data$get.sigma(pheno1, pheno2)

			status = list()
			output = data.frame(
				phenotype1 = pheno1, phenotype2 = pheno2,
				pair.type = private$data$get.types(pheno1, pheno2)[1,2],
				no.snps = private$data$no.snps(), no.components = no.components,
				omega.var1 = omega[1,1], omega.var2 = omega[2,2],	omega.cov = omega[1,2],
				rho = private$data$get.rg(pheno1, pheno2)[1,2], r2 = NA,
				stringsAsFactors=F
			)

			if (private$settings$get("compute.ci") != "never") output[c("rho.lower", "rho.upper", "r2.lower", "r2.upper")] = NA
			if (private$settings$get("compute.pval")) {
				output$test.type = ifelse(output$pair.type == "truncated", "normal", private$settings$get("use.method"))
				output[c(ifelse(output$test.type == "sampling", "no.iter", "sampling.variance"), "p.value")] = NA
			}

			if (!is.na(output$rho)) {
				if (abs(output$rho) > private$settings$get("rg.limit")) output$rho = NA
				output$r2 = output$rho^2
				if (is.na(output$rho)) status$estimation = "estimate out of bounds"
			} else status$model = ifelse(grepl("^unavailable", output$pair.type), gsub("unavailable - ", "", output$pair.type), "")

			if (!is.na(output$rho)) {
				if (private$settings$get("compute.pval")) {
					if (output$test.type == "sampling") {
						res = muffleErrors(integral.p(bivariate.integral, K=no.components, omega=omega/no.components, sigma=sigma, adap.thresh=private$settings$get("adaptive.threshold")))
						if (!is.error(res)) output[c("no.iter", "p.value")] = res[c("no.iter", "p.value")]
						else status$p.value = error.message(res)
					} else {
	# TODO: implement skew normal
						output$sampling.variance = private$data$sampling.variance(pheno1, pheno2)
						output$p.value = pnorm(-abs(omega[1,2] / sqrt(output$sampling.variance)))*2
					}
				}

				if (private$settings$get("compute.ci") == "always" || (private$settings$get("compute.ci") == "significant" && !is.na(output$p.value) && output$p.value <= private$settings$get("ci.threshold"))) {
					if (output$pair.type != "truncated") {
						res = muffleErrors(ci.bivariate(K=no.components, omega=omega/no.components, sigma=sigma))
						if (!is.error(res)) {
							output[c("rho.lower", "rho.upper")] = res$rho
							output[c("r2.lower", "r2.upper")] = res$r2
						} else status$ci = error.message(reS)
					} else status$ci = "not available for truncated phenotypes"
				}
			}
			return(list(output=output, status=status))
		}
	),
	public = list(
		compute = function(targets=NULL, phenotypes=NULL) {
			phenotypes = private$resolve.phenotypes(phenotypes)
			targets = if (!is.null(targets)) private$resolve.phenotypes(targets) else phenotypes

			results = list(); status = list()
			for (phs in make.pairs(targets, phenotypes, allow.equal=F, ignore.order=T, sort=T)) {
				label = paste0(phs, collapse=" ~ ")
				curr = private$compute.rg(phs[1], phs[2])
				results[[label]] = data.frame(model.name=label, curr$output, stringsAsFactors=F)
				if (length(curr$status) > 0) status[[label]] = curr$status
			}

			results = add.rownames(fill.rowmerge(results), names=NULL)
			return(list(results=results, status=status))
		}
	)
)



