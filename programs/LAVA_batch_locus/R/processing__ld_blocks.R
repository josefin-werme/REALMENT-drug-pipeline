#' @include global.R
#' @include input__ld_reference.R



## LD and decomposition for block of SNPs
LDblock = R6::R6Class("LDblock",
	private = list(
		settings = DEFAULT.SETTINGS("LDblock",
			prune.threshold = add.info(globals$pca.prune.threshold, type="numeric", range=c(0.1, 100)),  ## proportion of variance to retain after decomposition; can be on 0-1 or 0-100 scale
			max.components = add.info(Inf, type="numeric", min=1)   ## upper bound on principal components, excess is discarded
		),
		source.info = list(id=NULL, label=NULL),
		block.id = NULL,
		snp.index = NULL,  ## SNPindex object
		snp.info = NULL, ## data.frame with columns internal.id, snp.id, chromosome, position, allele1, allele2
		decomposition = list(no.components=NULL, vectors=NULL, values=NULL, components=NULL),  ## pruned decomposition of LD matrix

		get.printer = function() {
			printer = ObjectPrinter$new("LD Block")$parameter.list(source=private$source.info$label)
			no.comp = self$no.components()
			printer$add.parameter("no. variants", self$no.snps())$add.parameter("no. components", ifelse(is.na(no.comp), "unknown", no.comp))
			return(printer)
		},

		validate.settings = function(...) {
			private$settings = initialize.settings(self, ...)
			if (private$settings$get("prune.threshold") > 1) private$settings$set.values(prune.threshold = private$settings$get("prune.threshold")/100)
		},

		compute.decomposition = function() {undefined.error("compute.decomposition")},

		## subsets block to only those SNPs specified in internal.ids argument (discarding any not found in block)
		slice = function(internal.ids, decomposition=T) {undefined.error("slice", class(self)[1])}
	),
	public = list(
		initialize = function(source, ...) {
			private$validate.settings(...)
			private$source.info = list(id = source$data.id(), label = source$abbreviate())
		},
		print = function(...) {cat(private$get.printer()$to.string())},

		id = function() {if (is.null(private$block.id)) private$block.id = rlang::hash(c(private$source.info$id, private$snp.index$id(), self$get.decomposition()$values)); return(private$block.id)},
		equals = function(other) {return(class(self)[1] == class(other)[1] && self$id() == other$id())},

		no.snps = function() {return(private$snp.index$size())},
		no.components = function() {return(self$get.decomposition()$no.components)},

		get.snp.index = function() {return(private$snp.index)},
		get.snp.info = function(columns=NULL, collapse=T) {
			info = if (!is.null(columns)) private$snp.info[columns] else private$snp.info
			return(if (collapse && ncol(info) == 1) info[,1] else info)
		},

		get.decomposition = function() {if (is.null(private$decomposition$no.components)) private$compute.decomposition(); return(private$decomposition)},
		get.vectors = function() {return(self$get.decomposition()$vectors)}, # Q
		get.values = function() {return(self$get.decomposition()$values)}, # L
		get.root = function() {return(self$get.decomposition()$vectors %*% diag(sqrt(self$get.decomposition()$values)))}, ## Q %*% diag(sqrt(L))
		get.inverse.root = function() {return(self$get.decomposition()$vectors %*% diag(sqrt(1/self$get.decomposition()$values)))}, ## Q %*% diag(sqrt(1/L))

		## create LD block containing subset of SNPs based on internal.ids (either SNPindex object or vector of internal IDs)
		subset = function(internal.ids, decomposition=T) {
			if (!has.type(internal.ids, "SNPindex")) internal.ids = SNPindex$new(internal.ids)
			if (!private$snp.index$contains(internal.ids)) input.error("LD block does not contain all requested SNPs, cannot retrieve subblock")
			return(self$clone()$.__enclos_env__$private$slice(internal.ids, decomposition=decomposition))
		},

		## create LD block containing subset of SNPs shared with internal.ids argument
		intersection = function(internal.ids, decomposition=T) {return(self$clone()$.__enclos_env__$private$slice(internal.ids, decomposition=decomposition))},

		## check/align data in data.frame with internal.id column to LD block,
		## - mode=truncate checks that all SNPs are in data, discards excess in data and reorders if needed
		## - mode=sort checks that SNPs are the same, but reorders data if needed
		## - mode=strict checks that SNPs are the same and already in same order, error if not
		## if the function returns without error, the input data will be aligned to the LD block
		align.data = function(data, mode=c("truncate", "sort", "strict")) {mode = match.arg(mode)
			if (!is.data.frame(data) || !("internal.id" %in% names(data))) fatal.error("invalid input to align.data()")
			if (mode == "truncate") data = data[data$internal.id %in% private$snp.index$snps(),]
			if (mode %in% c("sort", "truncate")) data = data[order(data$internal.id),]
			if (!private$snp.index$equals(data$internal.id)) input.error("data is not aligned with LD block")
			return(data)
		}
	)
)

## LD block based on raw genotype input data
LDblockRaw = R6::R6Class("LDblockRaw",
	inherit = LDblock,
	private = list(
		data = NULL,  ## data.frame with standardized genotype data (missing values imputed to 0)
		indiv.info = NULL,  ## data.frame with indiv.id, family.id, sex, phenotype (if loaded)

		process.data = function(internal.ids, data) {
			if (ncol(data) != internal.ids$size()) input.error("size of genotype data does not match internal ID vector")
			if (!is.null(private$indiv.info) && nrow(data) != nrow(private$indiv.info)) input.error("size of genotype data does not match provided individual-level information")

			data = sweep(data.matrix(data), 2, apply(data, 2, mean, na.rm=T), FUN="-")
			data[is.na(data)] = 0
			sd = apply(data, 2, sd, na.rm=T)
			private$data = sweep(data, 2, sd, FUN="/")

			invalid = is.na(sd) | sd <= 0
			if (any(invalid)) {
				data = data[,!invalid, drop=F]
				internal.ids = SNPindex$new(internal.ids$snps[!invalid])
			}
			private$snp.index = internal.ids
		},

		compute.decomposition = function() {
			if (ncol(private$data) > 1) {
				decomp = try(svd(private$data), silent=T)
				if (class(decomp) == "try-error") decomp = try(decomp(private$data), silent=T)
				if (class(decomp) == "try-error") {
					decomp = eigen(cov(private$data))
					vectors = decomp$vectors
					values = decomp$values
				} else {
					vectors = decomp$v
					values = decomp$d * decomp$d / (nrow(private$data) - 1) * sign(decomp$d)
					values = values[decomp$d > 0]
				}

				req.comp = which(cumsum(values / sum(values)) >= private$settings$get("prune.threshold") & values > 0)
				no.components = min(req.comp, length(values), private$settings$get("max.components"))
				private$decomposition = list(no.components = no.components, vectors = vectors[,1:no.components], values = values[1:no.components])
			} else {
				if (ncol(private$data) == 1) private$decomposition = list(no.components = 1, vectors = matrix(1, 1, 1), values = 1)
				else private$decomposition = list(no.components = 0, vectors = matrix(0, 0, 0), values = numeric(0))
			}
		},

		## modifies object, should only be called directly after cloning
		## subsets block to only those SNPs specified in internal.ids argument (discarding any not found in block)
		slice = function(internal.ids, decomposition=T) {
			if (!has.type(internal.ids, "SNPindex")) internal.ids = SNPindex$new(internal.ids)
			if (!private$snp.index$equals(internal.ids)) {
				subset = internal.ids$extract.index(private$snp.index)
				private$snp.info = private$snp.info[subset,]
				private$snp.index = SNPindex$new(private$snp.index$snps()[subset])
				private$data = private$data[,subset,drop=F]
				private$decomposition = NULL; private$block.id = NULL
				if (decomposition) private$compute.decomposition()
			}
			return(self)
		},

		get.printer = function() {return(super$get.printer()$set.name("LD Block (raw data)"))}
	),
	public = list(
		initialize = function(source, internal.ids, data, ..., indiv.info=NULL, decompose=T) {
			super$initialize(source, ...)
			check.types(internal.ids="SNPindex")

			private$process.data(internal.ids, data)
			if (decompose) private$compute.decomposition()

			private$snp.info = subset(source$get.snp.info(private$snp.index), select=-internal.id)
			private$indiv.info = indiv.info
		},

		no.indiv = function() {return(nrow(private$data))},
		get.indiv.info = function() {return(private$indiv.info)},

		get.data = function() {return(private$data)},

		get.components = function() {
			if (is.null(private$decomposition$components)) private$decomposition$components = private$data %*% self$get.inverse.root()
			return(private$decomposition$components)
		},



#TODO add source ID check here; + move to external?
		## compute correlations of genetic components with those in other LDblock
		component.correlations = function(other) {
			if (other$no.indiv() != self$no.indiv()) input.error("LD blocks are incompatible, cannot compute correlations")
			return(t(self$get.components()) %*% other$get.components() / (self$no.indiv() - 1))
		},

		trim.components = function(drop=NULL, truncate=NULL) {
			curr.pcs = self$no.components()
			if (!is.null(drop)) {
				if (!is.null(truncate)) input.error("for function trim.components(), cannot use 'truncate' argument if 'drop' argument is specified")
				if (min(drop) <= 0 || max(drop) > curr.pcs) input.error("for function trim.components(), elements of 'drop' argument are out of bounds")
				drop = 1:curr.pcs %in% drop
			} else drop = 1:curr.pcs > truncate
			if (!is.null(drop)) {
				private$decomposition$vectors = private$decomposition$vectors[,!drop]
				private$decomposition$values = private$decomposition$values[!drop]
				private$decomposition$no.components = length(private$decomposition$values)
				if (!is.null(private$decomposition$components)) private$decomposition$components = private$decomposition$components[,!drop]
				private$block.id = NULL
			}
			invisible(self)
		}
	)
)
