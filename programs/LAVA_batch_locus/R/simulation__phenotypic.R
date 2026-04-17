#' @include global.R
#' @include simulation__genetic.R



## outputs marginal SNP test statistics for P phenotypes
SimulationPhenotypicModel = R6::R6Class("SimulationPhenotypicModel",
	private = list(
		genetic.model = NULL,
		parameters = list(
			sample.size = NULL,   ## length P numeric vector
			heritability = NULL,   ## length P numeric vector
			sampling.correlation = NULL   ## CovarianceMatrix object
		),
		current = list(),

		generate.internal = function() {undefined.error("generate", class(self)[1])}
	),
	public = list(
		initialize = function(genetic.model, sample.size, heritability, correlation=NULL) {
			check.types(genetic.model="SimulationGeneticModel"); check.types(correlation="CovarianceMatrix", .allow.null=T)
			private$genetic.model = genetic.model; no.pheno = genetic.model$no.pheno()

			if (length(sample.size) != no.pheno) input.error("length of 'sample.size' argument does not match number of phenotypes")
			if (any(sample.size < globals$simulation.minimum.N)) input.error("simulation sample sizes should all be at least ", globals$simulation.minimum.N)
			private$parameters$sample.size = sample.size

			if (length(heritability) != no.pheno) input.error("length of 'heritability' argument does not match number of phenotypes")
			if (any(heritability < 0 | heritability > 0.99)) input.error("simulation heritability values should all be between 0 and 0.99")
			private$parameters$heritability = heritability

			phenotypes = private$genetic.model$phenotypes()
			if (!is.null(correlation)) {
				IF.ANY(missing = !(phenotypes %in% correlation$get.names()), THEN=input.error("sampling correlation matrix does not contain %phenotype% ", items.and=phenotypes[missing]))
				private$parameters$sampling.correlation = correlation$subset(phenotypes)$standardize()
			} else private$parameters$sampling.correlation = as.covariance(diag(length(sample.size)), names=phenotypes)
		},
		clear = function() {private$current = list()},

		generate = function(redraw.genetic=T) {
			if (is.null(private$genetic.model)) input.error("cannot generate draw, no genetic model has been provided")
			if (redraw.genetic) private$genetic.model$generate()
			return(add.names(as.data.frame(private$generate.internal()), names=private$genetic.model$phenotypes()))
		},

		no.pheno = function() {return(length(private$parameters$sample.size))},

		get.genetic = function() {return(private$genetic.model$clone())},
		get.parameters = function(incl.genetic=T) {return(if (incl.genetic && !is.null(private$genetic.model)) list.merge(private$genetic.model$get.parameters(), private$parameters) else private$parameters)},

		variance = function() {undefined.error("variance", class(self)[1])},
		sigma = function() {
			variance = self$variance() * (1 - private$parameters$heritability) / (private$parameters$sample.size - 1)
			return(private$parameters$sampling.correlation$scale(variance))
		}
	)
)


## simulates normally distributed phenotypes with variance of 1
SimulationContinuousPhenotype = R6::R6Class("SimulationContinuousPhenotype",
	inherit = SimulationPhenotypicModel,
	private = list(
		generate.internal = function() {
			residual = private$current$ld.root %*% private$current$sampler$generate()
print(cor(residual))
			beta = sweep(private$genetic.model$beta(), 2, private$current$genetic.rescale, FUN="*") + residual
			samp.var = sweep(sweep(-beta^2, 2, self$variance(), FUN="+"), 2, (private$parameters$sample.size - 2), FUN="/")
			samp.var[samp.var < 0] = NA
			return(beta / sqrt(samp.var))
		}
	),
	public = list(
		initialize = function(genetic.model, sample.size, heritability, correlation=NULL) {
			super$initialize(genetic.model, sample.size, heritability, correlation)

			private$current$genetic.rescale = sqrt(private$parameters$heritability / private$genetic.model$variance())
			private$current$ld.root = private$genetic.model$get.ld()$get.root()
			private$current$sampler = MatrixSampler$new(self$sigma(), ncol(private$current$ld.root))
		},

		variance = function() {rep(1, self$no.pheno())}
	)
)



