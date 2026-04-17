#' @include global.R
#' @include processing__locus.R
#' @include model__sampling.R
#' @include model__confidence_intervals.R




## given a ProcessedLocus object, exposes cov(G) estimates and related moments, as well as distribution generators
## input ProcessedLocus object is stripped of all unusable phenotypes upon initialization
LocusModel = R6::R6Class("LocusModel",
	private = list(
		ld.block = NULL, ## global LD block for locus
		marginal = list(), ## named list of MarginalEstimates objects
		pair.type = matrix(), ## named symmetric matrix specifying type for each phenotype pair, currently: "harmonized", "truncated", "unavailable: [reason]"
		estimates = list(sigma = NULL, omega = NULL),
		cache = list(component.products=list()),

		pair.id = function(phenotypes) {return(paste0(sort(phenotypes), collapse=" - "))},

		## product = M = cov(W1,W2); trace = tr(t(M) %*% M)
		component.product = function(phenotypes, compute.trace=F) {
			marginal = private$marginal[phenotypes]
			id = private$pair.id(phenotypes)

			if (!is.null(private$cache$component.product[[id]])) {
				cross = private$cache$component.product[[id]]
				if (nrow(cross$product) != marginal[[1]]$no.components()) cross$product = t(cross$product)
			} else {
				cross = marginal[[1]]$get.ld()$component.correlations(marginal[[2]]$get.ld())
				private$cache$component.product[[id]] = list(product = cross)
			}

			if (compute.trace && is.null(cross$trace)) {
				cross$trace = sum(diag(cross$product %*% t(cross$product)))
				private$cache$component.product[[id]]$trace = cross$trace
			}

			return(cross)
		},

		## match partial names, check and throw error for ambiguous and unknown, remove duplicates/expand wildcards
		check.phenotypes = function(...) {phenotypes = c(...)
			if (length(phenotypes) > 0) return(check.phenotypes(phenotypes, names(private$marginal), label="locus model"))
			else if (is.null(phenotypes)) return(NULL)
			else return(character(0))
		},

		check.pair = function(ph1, ph2) {
			phenotypes = private$check.phenotypes(ph1, ph2)
			if (length(phenotypes) != 2) input.error(ifelse(length(phenotypes) > 2, "more", "fewer"), " than two phenotypes specified")
			return(phenotypes)
		},

		process = function(data) {
			private$ld.block = data$get.ld()
			private$marginal = data$get.estimates()

			if (length(private$marginal) > 1) {
				sigma.scale = sqrt(unlist(lapply(private$marginal, function(est) {est$get("sigma")})))
				sigma.corr = data$get.correlations()$get()

				private$estimates$sigma = sweep(sweep(sigma.corr, 1, sigma.scale, FUN="*"), 2, sigma.scale, FUN="*")
				private$estimates$omega = diag.matrix(sapply(private$marginal, function(est) {est$get("omega")}), add.names=names(private$marginal))

				no.pheno = length(private$marginal)
				private$pair.type = diag.matrix(rep(NA, no.pheno), add.names=names(private$marginal))
				if (no.pheno > 1) {
					for (i in 1:(no.pheno-1)) {
						for (j in (i+1):no.pheno) {
							res = private$process.pair(names(private$marginal)[c(i,j)])
							private$estimates$omega[i,j] = private$estimates$omega[j,i] = res$omega
							private$pair.type[i,j] = private$pair.type[j,i] = res$type
						}
					}
				}
			} else if (length(private$marginal) == 1) {
				private$estimates$sigma = add.names(matrix(private$marginal[[1]]$get("sigma")), names=names(private$marginal))
				private$estimates$omega = add.names(matrix(private$marginal[[1]]$get("omega")), names=names(private$marginal))
				private$pair.type = add.names(matrix(NA,1,1), names=names(private$marginal))
			} else private$estimates$sigma = private$estimates$omega = private$pair.type = matrix(0,0,0)
		},

		## determines pair type and computes genetic covariance for pair of input phenotypes
		process.pair = function(phenotypes) {
			marginal = private$marginal[phenotypes]
			sigma = private$estimates$sigma[phenotypes, phenotypes]
			delta = lapply(marginal, function(est) {est$get("delta")})
			K = sapply(marginal, function(est) {est$no.components()})

			if (is.na(sigma[1,2])) return(list(omega=NA, type="unavailable - sample overlap is unknown"))

			global.ld = sapply(marginal, function(est) {est$get.ld()$equals(private$ld.block)})
			if (all(global.ld) || marginal[[1]]$get.ld()$equals(marginal[[2]]$get.ld())) { ## same SNPs and PCs
				omega = sum(delta[[1]] * delta[[2]]) - K[1] * sigma[1,2]
				return(list(omega=omega, type="harmonized"))
			} else if (any(global.ld)) {
				if (sigma[1,2] == 0) {
					cross = private$component.product(phenotypes)
					omega = t(delta[[1]]) %*% cross %*% delta[[2]]
					return(list(omega=omega, type="truncated"))
				} else return(list(omega=NA, type="unavailable - truncated phenotype with non-zero sample overlap"))
			} else return(list(omega=NA, type="unavailable - both phenotypes are truncated"))
		},

		compute.variance = function(phenotypes, omega.null=0) {
			marginal = private$marginal[phenotypes]
			sigma = private$estimates$sigma[phenotypes, phenotypes]
			omega = private$estimates$omega[phenotypes, phenotypes]
			delta = lapply(marginal, function(est) {est$get("delta")})
			K = sapply(marginal, function(est) {est$no.components()})

			type = self$get.types(phenotypes)[1,2]
			if (type == "harmonized") {
				var.base = K[1]*sigma[1,1]*sigma[2,2] + sigma[1,1]*omega[2,2] + sigma[2,2]*omega[1,1]
				var.overlap = 2*omega.null*sigma[1,2] + K[1]*sigma[1,2]^2
				return(as.numeric(var.base + var.overlap))
			} else if (type == "truncated") {
				cross = private$component.product(phenotypes, compute.trace=T)
				var.proj = c(
					sum((t(cross$product) %*% delta[[1]])^2) - sigma[1,1] * cross$trace,
					sum((cross$product %*% delta[[2]])^2) - sigma[2,2] * cross$trace
				)
				var.proj[var.proj < 0] = 0
				var.base = cross$trace*sigma[1,1]*sigma[2,2] + sigma[1,1]*var.proj[2] + sigma[2,2]*var.proj[1]
				return(var.base)
			} else return(NA)
		}
	),
	public = list(
		initialize = function(data) {
			check.types(data="ProcessedLocus")
			private$process(data$filter("available"))
		},

		no.pheno = function() {return(length(private$marginal))},
		phenotypes = function() {return(if (self$no.pheno() > 0) names(private$marginal) else character(0))},

		no.snps = function() {return(private$ld.block$no.snps())},
		no.components = function() {return(private$ld.block$no.components())},
		get.ld = function() {return(private$ld.block)},

		get.types = function(...) {if (length(c(...)) > 0) {phenotypes = private$check.phenotypes(...); return(private$pair.type[phenotypes,phenotypes,drop=F])} else return(private$pair.type)},
		get.sigma = function(...) {if (length(c(...)) > 0) {phenotypes = private$check.phenotypes(...); return(private$estimates$sigma[phenotypes,phenotypes,drop=F])} else return(private$estimates$sigma)},
		get.omega = function(...) {if (length(c(...)) > 0) {phenotypes = private$check.phenotypes(...); return(private$estimates$omega[phenotypes,phenotypes,drop=F])} else return(private$estimates$omega)},

		get.rg = function(..., missing.bound=NULL, truncate=F) {
			omega = self$get.omega(...)
			if (nrow(omega) > 0) {
				rg = cov2cor(omega)
				if (!is.null(missing.bound)) rg[row(rg) != col(rg) & !is.na(rg) & abs(rg) > missing.bound] = NA
				if (truncate) {index = !is.na(rg) & abs(rg) > 1; rg[index] = sign(rg[index])}
				return(rg)
			} else return(omega)
		},

		sampling.variance = function(ph1, ph2=NULL, omega.null=0) {return(private$compute.variance(private$check.pair(ph1, ph2), omega.null))}
	)
)



## base class for different types of analysis model

AnalysisModel = R6::R6Class("AnalysisModel",
	private = list(
		settings = DEFAULT.SETTINGS("AnalysisModel",
			compute.pval = T,  ## compute p-values yes/no
			compute.ci = c("never", "significant", "always"),  ## if "significant", will compute confidence intervals only if corresponding p-value is below ci.threshold parameter (sets compute.pval=T)
			ci.threshold = add.info(globals$conf.interval.threshold, type="numeric", range=c(0,1)),  ## p-value threshold at which to compute confidence intervals
			adaptive.threshold = add.info(globals$adaptive.pval.threshold, vector=T, type="numeric", range=c(0,1)),  ## vector of p-value thresholds to increase
			rg.limit = add.info(globals$rg.invalid.bounds, type="numeric", min=1),   ## local rg absolute bound (considered invalid if exceeded)
			use.method = c("normal", "sampler")   ## preferred method to use, will use other if not available
		),
		data = NULL,   ## LocusModel object

		## check validity of phenotypes, match partial phenotype names against available phenotypes, and remove duplicates/expand wildcards
		## phenotypes=NULL input will be expanded to list of all phenotypes
		## phenotypes specified in exclude are removed from selection if present
		## if invert=T, complement of selection relative to all phenotypes is returned (exluded phenotypes are removed -after- inverting)
		resolve.phenotypes = function(phenotypes, exclude=NULL, invert=F) {
			reference = private$data$phenotypes()
			phenotypes = if (!is.null(phenotypes)) check.phenotypes(phenotypes, reference, label="analysis input") else reference
			exclude = if (!is.null(exclude)) check.phenotypes(exclude, reference, label="analysis input") else character(0)

			if (invert) return(reference[!(reference %in% c(phenotypes, exclude))])
			else return(phenotypes[!(phenotypes %in% exclude)])
		}
	),
	public = list(
		initialize = function(data, ...) {
			if (has.type(data, "ProcessedLocus")) private$data = LocusModel$new(data)
			else if (has.type(data, "LocusModel")) private$data = data
			else input.error("input to AnalysisModel object should be of type ProcessedLocus or LocusModel")

			private$settings = initialize.settings(self, ...)
			if (private$settings$get("compute.ci") == "significant") private$settings$set(compute.pval=T)
		},

		no.pheno = function() {private$data$no.pheno()},
		phenotypes = function() {private$data$phenotypes()},

		get.settings = function(as.list=F) {return(if (as.list) private$settings$get.all() else private$settings$clone())}
	)
)





