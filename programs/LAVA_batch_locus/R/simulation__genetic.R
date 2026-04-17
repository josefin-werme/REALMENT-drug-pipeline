#' @include global.R
#' @include processing__locus.R


## STATUS: volatile (current values), shallow
SimulationGeneticModel = R6::R6Class("SimulationGeneticModel",
	private = list(
		settings = NULL,
		ld.info = list(
			ld.block = NULL,  ## LDblock object
			root = NULL,
			inverse.root = NULL  ## stored transposed
		),
		parameters = list(
			no.pheno = NULL,
			omega = NULL,   ## CovarianceMatrix object
			structure = NULL,   ## StructuralModel object, may remain NULL
			distribution = list(type = NULL)
		),
		current = list(),   ## most recently generated values

		get.values = function(parameter) {
			if (is.null(private$current[[parameter]])) {
				if (parameter == "delta") self$generate()
				else if (parameter == "alpha") private$current$alpha = private$ld.info$inverse.root %*% private$get.values("delta")
				else if (parameter == "beta") private$current$beta = private$ld.info$root %*% private$get.values("delta")
			}
			return(private$current[[parameter]])
		},

		generate.delta = function() {undefined.error("generate.delta", class(self)[1])}
	),
	public = list(
		initialize = function(names, ld.block, genetic.structure, ...) {
			private$settings = initialize.settings(self, ...)
			check.types(ld.block="LDblock", genetic.structure=c("CovarianceMatrix", "StructuralModel"))
			private$ld.info = list(ld.block = ld.block, root = ld.block$get.root(), inverse.root = t(ld.block$get.inverse.root()))

			omega = if (has.type(genetic.structure, "StructuralModel")) genetic.structure$get.correlations() else genetic.structure

			if (length(names) == 0) input.error("cannot set SimulationGeneticModel object to zero phenotypes")
			IF.ANY(duplicate=duplicated(names), THEN=input.error("input contains duplicate phenotype %name% ", items.and=unique(names[duplicate])))
			IF.ANY(missing=!(names %in% omega$get.names()), THEN=input.error("genetic covariance matrix is missing %phenotype% ", items.and=names[missing]))

			private$parameters$no.pheno = length(names)
			private$parameters$omega = omega$subset(names)
			if (has.type(genetic.structure, "StructuralModel")) private$parameters$structure = genetic.structure$set.selection(names)

			if (!is.invertible(private$parameters$omega)) input.error("genetic covariance matrix is not invertible")
		},
		clear = function() {private$current = list()},

		generate = function() {self$clear(); private$current$delta = private$generate.delta(); invisible(self)},

		phenotypes = function() {return(private$parameters$omega$get.names())},
		no.pheno = function() {return(private$parameters$no.pheno)},
		no.snps = function() {return(private$ld.info$ld.block$no.snps())},
		no.components = function() {return(private$ld.info$ld.block$no.components())},

		get.ld = function() {return(private$ld.info$ld.block)},
		get.parameters = function() {return(list.merge(settings.genetic=private$settings$clone(), private$parameters))},

		omega = function(as.matrix=F) {return(if (as.matrix) private$parameters$omega$get() else private$parameters$omega)},
		variance = function() {return(add.names(diag(private$parameters$omega$get()), names=NULL))},

		alpha = function() {return(private$get.values("alpha"))},
		beta = function() {return(private$get.values("beta"))},
		delta = function() {return(private$get.values("delta"))}
	)
)


HomogeneousGeneticModel = R6::R6Class("HomogeneousGeneticModel",
	inherit = SimulationGeneticModel,
	private = list(
		sampler = NULL,

		generate.delta = function() {return(private$sampler$generate())}
	),
	public = list(
		initialize = function(names, ld.block, genetic.structure, ..., parameter = c("alpha", "delta")) {parameter = match.arg(parameter)
			super$initialize(names, ld.block, genetic.structure, ...)
			private$parameters$distribution = list(type="homogenous", parameter=parameter)
			private$sampler = MatrixSampler$new(private$parameters$omega, private$ld.info$ld.block$no.components(), weights=if (parameter == "alpha") sqrt(private$ld.info$ld.block$get.values()), exact=T, mode="product")
		}
	)
)


SparseGeneticModel = R6::R6Class("SparseGeneticModel",
	inherit = SimulationGeneticModel,
	private = list(
		settings = DEFAULT.SETTINGS("SparseGeneticModel",
			maximum.tries = add.info(5, min=1),   ## number of tries allowed, in case of inversion error during generation
			minimum.causal = add.info(5, min=1),
			causal.ratio = add.info(2, min=1),   ## number of causal SNPs per component is max(minimum.causal, [no. exogenous components] * causal.ratio)
			block.size.ratio = add.info(5, min=2),  ## minimum ratio of total to causal SNPs in blocks
			block.spacing = add.info(1, min=0)   ## mininum number of unused blocks between causal blocks
		),

		pattern.matrix = NULL,

		generate.exogenous = function() {
			distr = private$parameters$distribution
			P = distr$no.exogenous; K.causal = distr$causal.snps

			## select spaced random blocks by distributing the excess blocks over intervals between blocks (+ before first/after last blocks)
			spacing = head(rmultinom(1, distr$block.count - distr$block.minimum, rep(1, P+1)) + c(0, rep(private$settings$get("block.spacing"), P)), -1)
			block.index = sample(cumsum(spacing) + 1:P)
			snp.offset = sample(distr$offset.range+1, 1) - 1

			R = t(private$ld.info$root); delta = NULL
			for (p in 1:P) {
				causal = snp.offset + (block.index[p] - 1) * distr$block.size + sort(sample(distr$block.size, K.causal))
				if (p > 1) {
					D = R[,causal] %*% matrix(rnorm(K.causal * p), ncol=p)
					wts = solve(t(delta) %*% D[,-1]) %*% t(delta) %*% D[,1]
					delta = cbind(delta, D[,1] - D[,-1] %*% wts)
				} else delta = R[,causal] %*% rnorm(K.causal)
			}

			delta = sweep(delta, 2, sqrt(apply(delta^2, 2, sum)), FUN="/")
			return(delta[,sample(P,P)])
		},

		generate.delta = function() {
			for (i in 1:private$settings$get("maximum.tries")) {
				delta.base = muffleErrors(private$generate.exogenous())
				if (!is.error(delta.base)) break
			}
			if (is.error(delta.base)) data.error("reached maximum number of tries trying to generate genetic components")

			return(delta.base %*% private$pattern.matrix)
		}
	),
	public = list(
		initialize = function(names, ld.block, genetic.structure, ...) {
			check.types(genetic.structure="StructuralModel")
			super$initialize(names, ld.block, genetic.structure, ...)

			distribution = list(type = "sparse")
			distribution$no.exogenous = private$parameters$structure$no.pheno("exogenous")
			distribution$causal.snps = max(private$settings$get("minimum.causal"), private$settings$get("causal.ratio") * distribution$no.exogenous)
			distribution$block.size = as.numeric(distribution$causal.snps * private$settings$get("block.size.ratio"))
			distribution$block.count = floor(self$no.snps() / distribution$block.size)
			distribution$block.minimum = as.numeric((1 + private$settings$get("block.spacing")) * (distribution$no.exogenous - 1) + 1)
			distribution$offset.range = self$no.snps() %% distribution$block.size
			if (distribution$block.count < distribution$block.minimum) input.error("at current settings number of SNPs (", self$no.snps(), ") is insufficient to accommodate required ", items.count=distribution$no.exogenous, " exogenous %phenotype%")
			private$parameters$distribution = distribution

			model = private$parameters$structure$get.model()
			private$pattern.matrix = model$correlations$get()[model$exogenous, self$phenotypes(),drop=F]
		}
	)
)







