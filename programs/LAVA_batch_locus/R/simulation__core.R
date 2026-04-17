#' @include global.R
#' @include processing__ld_blocks.R
#' @include processing__locus.R
#' @include simulation__genetic.R
#' @include simulation__phenotypic.R





SimulationGenerator = R6::R6Class("SimulationGenerator",
	private = list(
		settings = DEFAULT.SETTINGS("SimulationGenerator",
			IMPORT.DEFAULTS("LDblock", "LocusData"),
			iteration.prefix = "iter",
			maximum.snps = add.info(Inf, type="integer", min=100)
		),
		ld.reference = NULL,   ## LDreference object
		input = list(),
		data = list(
			iteration = 0,   ## index of current iteration
			ld.block = NULL,  ## LDblock object
			genetic = NULL,   ## SimulationGeneticModel object
			phenotypic = NULL,   ## SimulationPhenotypicModel object
			template = NULL   ## SimulatedLocusData object
		),

		stage.index = function(stage=c("uninitialized", "ld.block", "genetic", "phenotypic", "ready")) {return(match(match.arg(stage), eval(formals(sys.function())$stage)))},
		check.stage = function(stage, exact=F) {diff = private$get.stage() - ifelse(is.character(stage), private$stage.index(stage), stage); return(diff == 0 || (!exact && diff > 0))},
		get.stage = function() {return(max(1, private$input$stage))},
		set.stage = function(stage) {private$input$stage = ifelse(is.character(stage), private$stage.index(stage), stage)},
		unset.stage = function(stage) {private$input$stage = max(ifelse(is.character(stage), private$stage.index(stage), stage) - 1, 1)},

		set.input = function(..., stage="uninitialized") {args = list(...)
			for (i in seq_along(args)) private$input[[names(args)[i]]] = args[[i]]
			private$unset.stage(stage)
		},

		set.data = function(required="ready") {
			if (!private$check.stage(required)) {
				private$data$iteration = 0
				for (stage in (private$get.stage()+1):private$stage.index(required)) {
					if (stage == private$stage.index("ld.block")) {
						if (is.null(private$input$locus)) input.error("no locus has been set")
						snp.index = private$input$locus$map.index(private$ld.reference$get.snp.info())
						ld.block = private$ld.reference$get.block(snp.index, private$settings, decompose=F)

						if (ld.block$no.snps() > private$settings$get("maximum.snps")) ld.block = ld.block$subset(head(ld.block$get.snp.index()$snps(), private$settings$get("maximum.snps")))
						if (ld.block$no.snps() < private$settings$get("minimum.snps")) input.error("fewer than ", private$settings$get("minimum.snps"), " SNPs in locus genotype data")
						if (ld.block$no.components() < private$settings$get("minimum.components")) input.error("fewer than ", private$settings$get("minimum.components"), " genetic principal components in locus genotype data")
						private$data$ld.block = ld.block
					} else if (stage == private$stage.index("genetic")) {
						if (is.null(private$input$genetic.structure)) input.error("genetic structure has not been specified")
						if (private$input$distribution == "sparse") private$data$genetic = SparseGeneticModel$new(private$input$phenotypes, private$data$ld.block, private$input$genetic.structure, private$input$genetic.settings)
						else private$data$genetic = HomogeneousGeneticModel$new(private$input$phenotypes, private$data$ld.block, private$input$genetic.structure, private$input$genetic.settings, parameter=private$input$distribution)
					} else if (stage == private$stage.index("phenotypic")) {
						no.pheno = private$data$genetic$no.pheno()
						h2 = if (length(private$input$heritability) == 1) rep(private$input$heritability, no.pheno) else private$input$heritability
						if (length(h2) != no.pheno) input.error("number of heritability values does not match number of phenotypes")
						N = if (length(private$input$sample.size) == 1) rep(private$input$sample.size, no.pheno) else private$input$sample.size
						if (length(N) != no.pheno) input.error("number of sample size values does not match number of phenotypes")
						if (!is.null(private$input$sampling.correlations)) {
							IF.ANY(missing = !(private$input$phenotypes %in% private$input$sampling.correlations$get.names()), THEN=input.error("sampling correlation matrix is missing %phenotype% ", items.and=private$input$phenotypes[missing]))
							private$input$sampling.correlations = private$input$sampling.correlations$subset(private$input$phenotypes)
						}
						private$data$phenotypic = SimulationContinuousPhenotype$new(private$data$genetic, N, h2, private$input$sampling.correlations)

						pheno.info = lapply(1:no.pheno, function(i) {SimulatedPhenotype$new(private$input$phenotypes[i], sample.size=N[i], heritability=h2[i])})
						private$data$template = SimulationLocusData$new(private$input$locus, private$data$ld.block, pheno.info, private$data$phenotypic$sigma(), private$settings)
					}
					private$set.stage(stage)
				}
			}
		}
	),
	public = list(
		initialize = function(ld.reference, ...) {
			check.types(ld.reference="LDreferenceData"); private$ld.reference = ld.reference
			private$settings = initialize.settings(self, ...)
		},
		clear = function() {private$input = list(); private$data = list()},

		no.pheno = function() {private$set.data("genetic"); return(private$data$genetic$no.pheno())},
		no.snps = function() {private$set.data("ld.block"); return(private$data$ld.block$no.snps())},
		no.components = function() {private$set.data("ld.block"); return(private$data$ld.block$no.components())},

		get.parameters = function() {private$set.data("phenotypic"); return(list.merge(locus = private$input$locus, private$data$phenotypic$get.parameters()))},
		get.ld = function() {private$set.data("ld.block"); return(private$data$ld.block)},

		change.settings = function(...) {private$settings$set(...); private$set.stage("uninitialized"); invisible(self)},
		reset.count = function() {private$data$iteration = 0; invisible(self)},

		## locus argument should either be an existing locus object, or a locus name
		## in the latter case, additional arguments should contain either contain chromosome/start/stop, or a snps argument with a list of SNPs
		set.locus = function(locus, ..., maximum.snps=NULL) {args = list(...)
			if (!is.null(maximum.snps)) private$settings$set(maximum.snps=maximum.snps)
			if (!has.type(locus, "LocusDefinition")) locus = create.locus(locus=locus, ...)
			private$set.input(locus=locus)
			invisible(self)
		},

		## genetic.structure can be either a CovarianceMatrix (rescaled to correlation) or a StructuralModel
		set.genetic = function(genetic.structure, ..., distribution=c("alpha", "delta", "sparse"), use=NULL) {distribution = match.arg(distribution)
			check.types(genetic.structure=c("CovarianceMatrix", "StructuralModel"))

			genetic.structure = genetic.structure$clone(); use = unique(use)
			available = if (has.type(genetic.structure, "StructuralModel")) genetic.structure$phenotypes("selection") else genetic.structure$get.names()
			if (!is.null(use)) {
				IF.ANY(missing=!(use %in% available), THEN=input.error("specified genetic structure does not contain required %phenotype% ", items.and=use[missing]))
				if (has.type(genetic.structure, "StructuralModel")) genetic.structure$set.selection(use)
				else genetic.structure = genetic.structure$subset(use)
			} else use = available
			if (has.type(genetic.structure, "CovarianceMatrix")) genetic.structure = as.covariance(genetic.structure, type="genetic correlation matrix")$standardize()

			private$set.input(phenotypes=use, genetic.structure=genetic.structure, distribution=distribution, genetic.settings=list(...), clear="genetic")
			invisible(self)
		},

		## heritability and sample.size can be single value or P-length vector; currently, only simulating continuous phenotype
		set.phenotypic = function(heritability=NULL, sample.size=NULL, correlations=NULL) {
			check.types(correlations="CovarianceMatrix", .allow.null=T)
			if (!is.null(correlations)) correlations = correlations$set.type("sampling correlation matrix")$standardize()
			private$set.input(heritability=heritability, sample.size=sample.size, sampling.correlations=correlations, clear="phenotypic")
			invisible(self)
		},

		## only redraws residuals if redraw.genetic=F; resets iteration count back to 1 if reset.count=T
		generate = function(redraw.genetic=T, reset.count=F) {
			if (!private$check.stage("ready")) private$set.data()
			stats = private$data$phenotypic$generate(redraw.genetic=redraw.genetic)
			private$data$iteration = ifelse(reset.count, 1, private$data$iteration + 1)
			data = private$data$template$update.id(paste0(private$settings$get("iteration.prefix"), private$data$iteration))
			return(data$set.statistics(stats))
		}
	)
)


