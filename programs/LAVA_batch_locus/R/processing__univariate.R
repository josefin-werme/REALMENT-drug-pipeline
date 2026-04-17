#' @include global.R
#' @include input__phenotypes.R
#' @include processing__ld_blocks.R


## SumstatBlock objects contain summary statistics and LD reference, computes relevant summary statistics
## UnivariateEstimates wraps around SumstatBlock and fits univariate models


##############################


## call compute.marginal() to obtain corresponding MarginalEstimates object

## STATUS: stable, deep (LDblock)
SumstatBlock = R6::R6Class("SumstatBlock",
	private = list(
		pheno.info = NULL, ## PhenotypeInfo object
		ld.block = NULL, ## LDblock object
		snp.data = NULL,	## data.frame with processed summary statistics, harmonized with ld.block
		computed = list(), ## list containing already computed internal statistics (to avoid recomputation on successive calls)

		fatal = function(...) {fatal.error(..., .label=list(phenotype=private$pheno.info$name()))},

		## return values for parameters that can exist at both the global and per-SNP level
		## - automatic: return per-SNP vector if available, single global value otherwise
		## - global: always a single value, will compute mean of per-SNP values if available
		## - snp: always a vector of length no.snps, will repeat global value if no per-SNP level data available
		resolve.level = function(name, output.type=c("automatic", "global", "snp")) {
			output.type = match.arg(output.type)
			if (name %in% names(private$snp.data)) {
				if (output.type == "global") {
					if (!(name %in% private$computed)) private$computed[[name]] = mean(private$snp.data[[name]], na.rm=T)
					return(private$computed[[name]])
				} else return(private$snp.data[[name]])
			} else {
				if (output.type == "snp") return(rep(private$pheno.info$get.global(name), nrow(private$snp.data)))
				else return(private$pheno.info$get.global(name))
			}
		},

		compute.xty = function() {undefined.error("compute.xty", class(self)[1])},
		compute.covar = function() {undefined.error("compute.covar", class(self)[1])}
	),
	public = list(
		initialize = function(ld.block, snp.data, pheno.info) {
			private$ld.block = ld.block; private$pheno.info = pheno.info
			private$snp.data = tryCatch(ld.block$align.data(snp.data, mode="sort"), error=function(err) {private$fatal("invalid input to SumstatBlock object, ", err)})
		},

		print = function(...) {
			printer = ObjectPrinter$new(class(self)[1])
			printer$parameter.list(phenotype=self$phenotype())$add.parameter("input type", self$data.type(as.string=T))$add.parameter("no. snps", self$no.snps())
			cat(printer$to.string())
		},
		clear = function() {private$computed = list()},

		phenotype = function() {return(private$pheno.info$name())},
		data.type = function(as.string=F, sep=" ") {return(private$pheno.info$get.traits(as.string=as.string, sep=sep))},

		is = function(...) {return(private$pheno.info$is(...))},  ## query phenotype traits, eg. is("binary")
		phenotype.info = function() {return(private$pheno.info)},

		no.snps = function() {return(nrow(private$snp.data))},
		no.components = function() {return(private$ld.block$no.components())},
		pheno.variance = function() {undefined.error("pheno.variance", class(self)[1])},

		snp.id = function() {return(private$ld.block$get.snp.info("snp.id"))},
		internal.id = function() {return(private$snp.data$internal.id)},
		alleles = function() {return(private$ld.block$get.snp.info(c("allele1", "allele2")))},
		sample.size = function(output.type=c("automatic", "global", "snp")) {return(private$resolve.level("sample.size", output.type))},
		statistic = function() {return(private$snp.data$statistic)},
		p.value = function() {return(2*pnorm(abs(private$snp.data$statistic), lower.tail=F))},
		xty = function() {
			if (!("xty" %in% names(private$computed))) private$compute.xty()
			return(private$computed$xty)
		},
		covariance = function(standardize=T) {
			if (!("covar.raw" %in% names(private$computed))) private$compute.covar()
			covar = private$computed$covar.raw
			if (standardize && self$pheno.variance() != 1) covar = covar / sqrt(self$pheno.variance())
			return(covar)
		},

		get.ld = function() {return(private$ld.block)},
		get.snp.data = function() {return(cbind(private$ld.block$get.snp.info(c("snp.id", "allele1", "allele2")), subset(private$snp.data, select=-internal.id)))},
		get.pheno.info = function() {return(private$pheno.info)},

		compute.marginal = function() {undefined.error("compute.marginal", class(self)[1])}
	)
)

ContinuousSumstatBlock = R6::R6Class("ContinuousSumstatBlock",
	inherit = SumstatBlock,
	private = list(
		compute.xty = function() {private$computed$xty = self$covariance(standardize=F) * (self$sample.size() - 1)},
		compute.covar = function() {
			beta = private$snp.data$statistic / sqrt(private$snp.data$statistic^2 + self$sample.size() - 2)
			if (self$pheno.variance() != 1) beta = beta * sqrt(self$pheno.variance())
			private$computed$covar.raw = beta
		}
	),
	public = list(
		initialize = function(ld.block, snp.data, pheno.info) {
			if (pheno.info$is("binary")) private$fatal("binary phenotype input into ContinuousSumstatBlock object")
			super$initialize(ld.block, snp.data, pheno.info)
		},

		pheno.variance = function() {return(1)},
		compute.marginal = function() {return(ContinuousMarginal$new(self))}
	)
)


## currently requires raw genotype data from reference to run
BinarySumstatBlock = R6::R6Class("BinarySumstatBlock",
	inherit = SumstatBlock,
	private = list(
		compute.xty = function() {
			genotypes = private$ld.block$get.data()
			sample.size = self$sample.size(output.type="snp"); case.proportion = self$case.proportion(output.type="snp")

			private$computed$xty = rep(0, nrow(private$snp.data))
			for (i in 1:nrow(private$snp.data)) {
				beta = find.beta(genotypes[,i], case.proportion[i]*sample.size[i], sample.size[i], private$snp.data$statistic[i])
				mu = 1 / (1+exp(-(beta[1] + beta[2]*genotypes[,i])))
				private$computed$xty[i] = sum(genotypes[,i] * mu)
			}
		},
		compute.covar = function() {private$computed$covar.raw = self$get.xty() / (self$sample.size() - 1)}
	),
	public = list(
		initialize = function(ld.block, snp.data, pheno.info) {
			if (pheno.info$is("continuous")) private$fatal("non-binary phenotype input into BinarySumstatBlock object")
			super$initialize(ld.block, snp.data, pheno.info)
		},

		pheno.variance = function() {
			sample.size = self$sample.size(output.type="global"); proportion = self$case.proportion(output.type="global")
			return(proportion * (1 - proportion) * sample.size / (sample.size - 1))
		},
		prevalence = function() {return(private$pheno.info$get.global("prevalence"))},
		case.proportion = function(output.type=c("automatic", "global", "snp")) {return(private$resolve.level("case.proportion", output.type))},

		compute.marginal = function() {return(LogisticBinaryMarginal$new(self))}
	)
)



##############################


## estimates are scaled to a phenotypic variance of one for continuous and binary phenotypes
## for meta-analysed phenotypes, scale is maintained at one and h2 is provided (equal to h2.composite)
## for other composite phenotypes, scale is undefined, and only h2.composite (weighted mean of input h2 values) is provided, and only if weights are not of mixed sign


## STATUS: stable, deep (SumstatBlock)
MarginalEstimates = R6::R6Class("MarginalEstimates",
	private = list(
		pheno.info = NULL,   ## PhenotypeInfo object
		estimates = list(delta = NA, eta = NA, sigma = NA, omega = NA, h2 = NA, h2.latent = NA, h2.composite = NA, p.value = NA),
		traits = list(continuous=F, binary=F, composite=F, meta=F),

		failure = function(...) {data.error(..., .label=list(phenotype=private$info$phenotype))},

		compute.chisq = function(delta, sigma) {return(pchisq(ifelse(sigma > 0, sum(delta^2) / sigma, NA), df=length(delta), lower.tail=F))},
		compute.F = function(delta, sigma, sample.size, no.components) {
			if (sample.size - no.components - 1 <= 0) input.error("degrees of freedom for F-test are zero or negative")
			return(pf(ifelse(sigma > 0, sum(delta^2) / sigma / no.components, NA), no.components, sample.size - no.components - 1, lower.tail=F))
		},

		## all estimates should be set via this function to validate, and should be set at the same time (old values are reset to NA if called again)
		## automatically sets omega and p.value (using chi-square test) if not provided, and truncates h2 estimats to [0,1] range
		set.estimates = function(delta, sigma, ...) {
			private$estimates[names(private$estimates)] = NA
			private$estimates$delta = delta

			if (!is.na(sigma) && sigma > 0) {
				private$estimates$sigma = sigma
				private$estimates = list.merge(private$estimates, flatten.arglist(...))

				if (is.na(private$estimates$omega)) private$estimates$omega = sum(delta^2) - length(delta)*sigma
				if (is.na(private$estimates$p.value)) private$estimates$p.value = private$compute.chisq(delta, sigma)

				for (h2 in grep("^h2", names(private$estimates), value=T)) private$estimates[[h2]] = min(max(private$estimates[[h2]], 0), 1)
			}
		}
	),
	public = list(
		phenotype = function() {return(private$pheno.info$name())},
		phenotype.info = function() {return(private$pheno.info)},

		no.snps = function() {return(self$get.ld()$no.snps())},
		no.components = function() {return(self$get.ld()$no.components())},

		## for parameters: delta, sigma, eta, omega, h2, h2.latent, p.value; returns NA if not present
		get = function(parameter) {
			if (!(parameter %in% names(private$estimates))) fatal.error("trying to query unknown marginal estimate parameter '", parameter, "'")
			return(if (!is.null(private$estimates[[parameter]])) private$estimates[[parameter]] else NA)
		},
		get.estimates = function() {return(private$estimates)},
		get.ld = function() {undefined.error("get.ld", class(self)[1])},

		is = function(trait) {
			if (!(trait %in% names(private$traits))) fatal.error("trying to query unknown marginal estimate trait '", trait, "'")
			return(private$traits[[trait]])
		}
	)
)


MarginalDataEstimates = R6::R6Class("MarginalDataEstimates",
	inherit = MarginalEstimates,
	private = list(
		sum.stats = NULL, ## SumstatBlock object
		sample.size = NULL,

		estimate = function() {undefined.error("estimate", class(self)[1])}
	),
	public = list(
		initialize = function(sum.stats) {
			if ("BinaryMarginal" %in% class(self)) private$traits$binary = TRUE
			else private$traits$continuous = TRUE

			check.types(sum.stats="SumstatBlock")
			private$pheno.info = sum.stats$phenotype.info()
			private$sum.stats = sum.stats
			private$sample.size = sum.stats$sample.size(output.type="global")

			private$estimate()

			if (self$no.components() / private$sample.size > globals$h2.component.ratio) private$estimates[grep("^h2", names(private$estimates), value=T)] = NA
		},

		get.ld = function() {return(private$sum.stats$get.ld())}
	)
)


ContinuousMarginal = R6::R6Class("ContinuousMarginal",
	inherit = MarginalDataEstimates,
	private = list(
		estimate = function() {
			N = private$sample.size; K = self$no.components()
			delta = t(private$sum.stats$get.ld()$get.inverse.root()) %*% private$sum.stats$covariance(standardize=T) / sqrt(private$sum.stats$pheno.variance())
			eta = (1 - sum(delta^2)) * (N - 1) / (N - K - 1)
			sigma = eta / (N - 1)
			h2 = 1 - eta
			p.value = private$compute.F(delta, sigma, N, K)

			private$set.estimates(delta, sigma, eta=eta, h2=h2, p.value=p.value)
		}
	)
)


BinaryMarginal = R6::R6Class("BinaryMarginal", inherit = MarginalDataEstimates)

LogisticBinaryMarginal = R6::R6Class("LogisticBinaryMarginal",
	inherit = BinaryMarginal,
	private = list(
		estimate = function() {no.iter=25; svar.thresh=1.75 #leaving in hard-coded parameters for now, remove later
			ld.block = private$sum.stats$get.ld(); genotypes = ld.block$get.data(); components = ld.block$get.components()
			N = private$sample.size; N.ref = nrow(genotypes); N.case = N * private$sum.stats$case.proportion(output.type="global")
			K = self$no.components()

			wty0 = N.case * N.ref / N
			wty1 = t(ld.block$get.inverse.root()) %*% private$sum.stats$xty()

			comp.use = 1:K
			while (length(comp.use) > 1) {
				W = cbind(1, components[,comp.use])
				wty = c(wty0, wty1[comp.use])

				beta = c(log(N.case/(N - N.case)), rep(0, length(comp.use)))
				for (i in 1:no.iter) {
					mu = 1 / (1+exp(-W%*%beta))
					s = as.numeric(mu * (1-mu))
					wsw.inv = try(solve(t(W) %*% diag(s) %*% W), silent=T)
					if (class(wsw.inv)[1] == "try-error") private$failure("inversion error when fitting multiple logistic regression model")
					beta = beta + wsw.inv %*% (wty - t(W) %*% mu)
				}
				V = diag(wsw.inv)[-1]
				svar.ratio = max(V) / quantile(V, 0.5)	# using ratio of max to median; will normally be very close to 1
				if (svar.ratio < svar.thresh) break 	# if above threshold, drop the PC with highest sampling variance and rerun
				comp.use = comp.use[-which.max(V)]
			}
			if (length(comp.use) < K) private$failure("unstable components in multiple logistic regression model", .sub.type="dropped.components", .data=list(dropped=which(!(1:K %in% comp.use))))

			linear.variance = sum(wty1^2) / ((N - 1) * N.ref/N)^2   ## ie. var(G) under linear regression of Y on W, G = W*delta.linear
			pheno.variance = private$sum.stats$pheno.variance()

			delta = beta[-1,] / sqrt(pheno.variance)
			sigma = mean(V) * N.ref / N / pheno.variance
			h2 = 1 - (1 - linear.variance / pheno.variance) * (N - 1) / (N - K - 1)   ## observed h2

			prevalence = private$sum.stats$prevalence()
			h2.latent = if (!is.na(prevalence)) h2 * prevalence/(1-prevalence) / dnorm(qnorm(prevalence))^2

			private$set.estimates(delta, sigma, h2=h2, h2.latent=h2.latent)
		}
	)
)



## h2.composite is set to weighted mean of input h2's unless weights have mixed signs
## meta-analysis subtype sets h2 to h2.composite, h2 is NA in general case
## correlation matrix should contain other phenotypes to compute correlations with
## - correlations among input phenotypes cannot be NA
CompositeEstimates = R6::R6Class("CompositeEstimates",
	inherit = MarginalEstimates,
	private = list(
		marginal = list(),   ## named list of MarginalEstimates objects
		weights = NULL,  ## vector of numeric weights, standardized to sum(abs(weights)) = 1
		correlations = NULL,   ## CovarianceMatrix object of sampling correlations with other phenotypes

		estimate = function(correlations) {
			delta = do.call(cbind, lapply(private$marginal, function(est) {est$get("delta")})) %*% private$weights

			sigma.input = sapply(private$marginal, function(est) {est$get("sigma")})
			if (any(is.na(sigma.input) | sigma.input <= 0)) input.error("cannot create composite phenotypes, input contains negative or missing sampling variance values")
			corr.input = correlations$subset(names(private$marginal))
			sigma = t(private$weights) %*% corr.input$scale(sigma.input)$get() %*% private$weights

			if (all(private$weights >= 0) || all(private$weights <= 0)) {
				h2.parts = sapply(private$marginal, function(est) {est$get(ifelse(est$is("composite"), "h2.composite", "h2"))})
				h2 = list(h2.composite = sum(abs(private$weights) * h2.parts))
				if (self$is("meta")) h2$h2 = h2$h2.composite
			} else h2 = list()

			private$set.estimates(delta, sigma[1], h2)

			corr.values = correlations$expand(self$phenotype(), variances=1)$get()
			var = t(private$weights) %*% corr.input$get() %*% private$weights
			covar = correlations$get()[,names(private$marginal)] %*% private$weights
			corr.values[correlations$get.names(),self$phenotype()] = corr.values[self$phenotype(),correlations$get.names()] = covar / sqrt(diag(correlations$get()) * var[1])

			private$correlations = correlations$set.values(corr.values)
		}
	),
	public = list(
		initialize = function(name, marginal, correlations, weights=NULL) {
			private$traits$composite = TRUE

			if (length(marginal) == 0 || !is.list(marginal) || !all(sapply(marginal, has.type, "MarginalEstimates"))) input.error("input to CompositeEstimates object should be a list of MarginalDataEstimates objects")
			if (any(sapply(marginal, has.type, "CompositeEstimates") & !sapply(marginal, has.type, "MetaAnalysisEstimates"))) input.error("composite phenotype cannot contain other composite phenotypes, unless meta-analysed")
			if (length(marginal) > 1 && !all(sapply(marginal[-1], function(est) {est$get.ld()$equals(marginal[[1]]$get.ld())}))) input.error("cannot create composite phenotype for phenotypes aligned to different LD blocks")
			private$marginal = add.names(marginal, names=sapply(marginal, function(est) {est$phenotype()}))

			if (!all(names(private$marginal) %in% correlations$get.names())) input.error("phenotypes in sampling correlation matrix input for composite phenotypes is incomplete")
			if (any(is.na(correlations$subset(names(private$marginal), ignore.missing=T)$get()))) input.error("sampling correlation matrix input for composite phenotype contains missing values")

			if (length(weights) <= 1) weights = rep(ifelse(length(weights) == 1, weights, 1), length(private$marginal))
			if (length(weights) != length(private$marginal)) input.error("incorrect number of weights for composite phenotype")

			weights[is.na(weights) | !is.finite(weights)] = 0; weights = round(weights, digits=10)
			if (length(marginal) == 1 || sum(abs(weights) >= 0) <= 1) input.error("composite phenotype requires at least two input phenotypes with non-zero weights")

			drop = weights == 0
			private$weights = weights[!drop] / sum(abs(weights))
			if (any(drop)) private$marginal = private$marginal[!drop]

			private$pheno.info = CompositePhenotype$new(name, lapply(private$marginal, function(est) {est$phenotype.info()}), private$weights, is.meta=self$is("meta"))
			private$estimate(correlations)
		},

		phenotype.info = function() {return(private$pheno.info)},
		get.ld = function() {return(private$marginal[[1]]$get.ld())},

		get.correlations = function() {return(private$correlations)}
	)
)


MetaAnalysisEstimates = R6::R6Class("MetaAnalysisEstimates",
	inherit = CompositeEstimates,
	public = list(
		initialize = function(name, marginal, correlations) {
			private$traits$meta = TRUE
			sigma = if (length(marginal) > 0 && is.list(marginal)) unlist(lapply(marginal, function(est) {return(ifelse(has.type(est, "MarginalEstimates"), est$get("sigma"), NA))}))
			super$initialize(name, marginal, correlations, weights=1/sigma)
		}
	)
)


