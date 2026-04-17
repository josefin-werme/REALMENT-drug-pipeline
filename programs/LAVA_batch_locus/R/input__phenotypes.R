#' @include global.R

## helper class to bundle and validate phenotype info; external input needs to contain at least 'phenotype' and 'filename' arguments
## primary RawPhenotype class type reads in input file header and validates, and provides interface for reading in full file
PhenotypeInfo = R6::R6Class("PhenotypeInfo",
	private = list(
		pheno.name = NULL, raw.info = list(),
		traits = list(continuous=F, binary=F, composite=F, meta=F, multi=F, localized=F, simulated=F),  ## for globals$phenotype.properties added in constructor

		error = function(...) {input.error(..., .label=list(phenotype=self$name()))},

		get.value = function(parameter) {
			if (!(parameter %in% names(private$raw.info))) private$error("missing parameter '", parameter, "' for ", class(self)[1], " object")
			return(private$raw.info[[parameter]])
		},

		get.printer = function() {undefined.error("get.printer", class(self)[1])}
	),
	public = list(
		initialize = function(...) {
			private$raw.info = flatten.arglist(..., filter.NA=T)
			private$pheno.name = private$get.value("phenotype")
			private$traits[globals$phenotype.properties] = FALSE

			IF.ANY(unknown=!(names(private$raw.info) %in% names(globals$info.header.index)), THEN=private$error("unknown %parameter% ", items.and=names(private$raw.info)[unknown], " for ", class(self)[1], " object"))
		},

		abbreviate = function() {return(paste0(self$name(), "  [", self$get.traits(), "]"))},
		print = function(...) {cat(private$get.printer()$to.string())},

		name = function(base.only=F) {return(private$pheno.name)},
		get.traits = function(as.string=T, sep=" ") {
			traits = names(private$traits)[unlist(private$traits)]
			if (as.string) return(ifelse(length(traits) > 0, paste(traits, collapse=sep), "generic"))
			else return(traits)
		},

		is = function(...) {traits = c(...)
			IF.ANY(unknown=!(traits %in% names(private$traits)), THEN=fatal.error("trying to query unknown phenotype %trait% ", items.and=traits[unknown]))
			return(unlist(private$traits[traits]))
		},

		get.global = function(name) {return(NA)}
	)
)


## info for phenotype directly read from data
## validation of input parameters is performed *here* (SummaryStatistics objects assume them to be valid without further checks), with get.configuration() exposing a list with all available parameters
## underlying data is exposed via get.data()

## possible configuration components
## - core: snp.id, allele1, allele2
## - metrics: sample.size, no.cases, no.controls, case.proportion (max. 2)
## - statistics: statistic -or- p.value, beta/log.odds/odds.ratio
## - localized: [localized phenotype units]

## valid binary column configurations (uses first available option in order listed here; other columns and global values are ignored)
## - two columns (will compute and store sample.size and case.proportion):
##   - sample.size and case.proportion
##   - no.cases and no.controls
##   - sample.size and one of no.cases/no.controls
## - sample size column plus global value(s)
##   - global case.proportion
##   - global no.cases and no.controls (will store computed case.proportion)
## - two globals (will compute and store sample.size and case.proportion):
##   - sample.size and case.proportion
##   - no.cases and no.controls
##   - sample.size and one of no.cases/no.controls
RawPhenotype = R6::R6Class("RawPhenotype",
	inherit = PhenotypeInfo,
	private = list(
		input.file = NULL,
		meta.data = list(
			sample.metrics = list(sample.size=NULL, case.proportion=NULL, prevalence=NULL),  ## global values, if provided and not available through input column in data
			localized = list(names=NULL, columns=list(), values=list()),    ## units for localized/multi phenotypes
			multi = list(name=NULL, columns=list(), value=NULL)  ## columns is named list, with vector of column names per sub-phenotype
		),

		configuration = list(core=NULL, metrics=NULL, statistics=NULL, units=NULL), ## list with available parameters in input data, by category (see comment above for possible parameters)
		data.interface = list(main=NULL, chromosomes=list()),  ## TextFileInputInterface object for input file, separate objects per chromosome if split

		## define additional parameters for InputInterface objects
		add.parameters = function(update.params) {
			if (length(update.params) > 0) {
				private$data.interface$main$update.parameters(update.params, append.existing=T)
				for (i in seq_along(private$data.interface$chromosomes)) {
					if (!is.null(private$data.interface$chromosomes[[i]])) private$data.interface$chromosomes[[i]]$update.parameters(update.params, append.existing=T)
				}
			}
		},

		## validate input file/files (expanding [CHR] placeholder) and create TextFileInputInterface object(s)
		## TextFileInputInterface objects read in file header for validation of required columns, and provides interface to full data file object
		validate.files = function() {
			private$input.file = private$get.value("filename")

			check.files.exist(private$input.file, resolve.chr=T)
			if (self$by.chromosome()) {
				chr.files = chromosome.files(private$input.file, prune.missing=T); columns = list()
				for (i in 1:nrow(chr.files)) {
					curr = TextFileInputInterface$new(chr.files$file[i], globals$sumstats.header.index)
					private$data.interface$chromosomes[[chr.files$chromosome[i]]] = curr
					columns[[i]] = sort(unique(tolower(curr$get.variables())))
				}
				if (length(columns) > 1 && !all(sapply(columns[-1], function(col) {length(col) == length(columns[[1]]) && all(col == columns[[1]])}))) private$error("headers of chromosome files for ", private$input.file, " are not consistent")
				private$data.interface$main = private$data.interface$chromosomes[[chr.files$chromosome[1]]]  ## set to first chromosome for subsequent column validation
			} else private$data.interface$main = TextFileInputInterface$new(private$input.file, globals$sumstats.header.index)
		},

		## load in and check general properties ('properties' parameter)
		validate.properties = function() {
			recognized = globals$phenotype.properties
			properties = tryCatch(parse.value.string(private$raw.info$properties), error=function(err) {private$error(err)})
			IF.ANY(unknown=!(names(properties) %in% recognized), THEN=private$error("%property% ", items.and=names(properties)[unknown], " specified in 'properties' %is% unknown"))

			for (prop in names(properties)) {if (prop %in% names(private$traits)) private$traits[[prop]] = TRUE}
			if ("truncated" %in% names(properties) && !check.mode("beta")) fatal.error("support for truncated phenotypes is currently disabled")
			if ("multi" %in% names(properties)) private$meta.data$multi$name = properties$multi
		},

		## validate user-defined input columns
		validate.columns = function() {
			update.params = tryCatch(parse.value.string(private$raw.info$parameters, allow.unnamed=F, infer.names=F), error=function(err) {private$error(err)})
			IF.ANY(unknown=!private$data.interface$main$has.parameters(names(update.params), aggregate=F), THEN=private$error("%parameter% ", items.and=names(update.params)[unknown], " specified in 'parameters' %is% unknown"))
			private$add.parameters(update.params)
		},

		## validate localized phenotype unit settings, and define corresponding parameters in InputInterface
		validate.localized = function() {
			private$traits$localized = TRUE
			units = tryCatch(parse.value.string(private$raw.info$unit, allow.unnamed=F, infer.names=T), error=function(err) {private$error(err)})
			units = lapply(units, function(u) {unique(tolower(u))})

			if (self$is("multi") && any(names(units) == private$meta.data$multi$name)) private$error("unit name '", private$meta.data$multi$name, "' conflicts with multi-input type name")
			IF.ANY(duplicate.names=duplicated(names(units)), THEN=private$error("unit %name% ", items.and=unique(names(units)[duplicate.names]), " %is% used more than once"))
			IF.ANY(duplicate.columns=duplicated(unlist(units)), THEN=private$error("%column% ", items.and=unique(unlist(units)[duplicate.columns]), " %is% %[+each]% specified for multiple units"))

			IF.ANY(exist=private$data.interface$main$update.parameters(units, append.existing=F), THEN=private$error("parameter %name% ", items.and=exist, " %is% already in use"))
			IF.ANY(unavailable=!private$data.interface$main$has.available(names(units)), THEN=private$error("%column% specified for %unit% ", items.and=names(units)[unavailable], " %does% not exist"))

			private$add.parameters(units)
			private$meta.data$localized = list(names=names(units), columns=units, values=list())
			private$configuration$localized = names(units)
		},

		core.columns = function() {
			## validate core parameters
			private$configuration$core = c("snp.id", if (!self$is("prealigned")) c("allele1", "allele2"))
			IF.ANY(unavailable=!private$data.interface$main$has.available(private$configuration$core), THEN=private$error("%column% for required %parameter% ", items.and=private$configuration$core[unavailable], " %is% missing"))
		},

		statistic.columns = function() {
			effect.params = c(if (self$is("binary")) c("log.odds", "odds.ratio"), "beta")
			if (self$is("multi")) {
				columns = list(); add.index = list()
				available = private$data.interface$main$map.partial(c("statistic", "p.value", effect.params))
				available = available[order(match(available$parameter, c("statistic", "p.value", effect.params))),]
				available$new.param = paste0(available$parameter, "::", available$label)
				for (label in unique(available$label)) {
					curr = available[available$label == label,]
					if (curr$parameter[1] == "statistic" || (curr$parameter[1] == "p.value" && nrow(curr) >= 2)) {
						use = if (curr$parameter[1] == "p.value") 1:2 else 1
						columns[[label]] = curr$new.param[use]
						add.index = c(add.index, named.list(names=curr$new.param[use], curr$input.name[use]))
					}
				}
				if (length(columns) == 0) private$error("unable to identify phenotypes in multi-phenotype input file")

				private$configuration$statistics = names(add.index)
				private$meta.data$multi$columns = columns
				private$add.parameters(add.index)
			} else {
				if (!private$data.interface$main$has.available("statistic")) {
					if (private$data.interface$main$has.available("p.value")) {
						available = private$data.interface$main$has.available(effect.params)
						if (any(available)) private$configuration$statistics = c("p.value", effect.params[available][1])
						else private$error("column for %[+one of]% %parameter% ", items.or=effect.params, " is required if no column for parameter 'statistic' is present")
					} else private$error("column for parameter 'p.value' is required if no column for parameter 'statistic' is present")
				} else private$configuration$statistics = "statistic"
			}
		},

		metric.columns = function() {
			if (self$is("binary")) {
				binary.params = c("sample.size", "no.cases", "no.controls", "case.proportion")
				has.column = add.names(private$data.interface$main$has.available(binary.params), binary.params)

				if (has.column["sample.size"] && has.column["case.proportion"]) private$configuration$metrics = c("sample.size", "case.proportion")
				else if (has.column["no.cases"] && has.column["no.controls"]) private$configuration$metrics = c("no.cases", "no.controls")
				else if (has.column["sample.size"] && any(has.column[2:3])) private$configuration$metrics = head(binary.params[has.column], 2)   ## use sample.size and one of no.cases or no.controls
				else {  ## don't have enough information in input file columns, add global values
					has.global = binary.params %in% names(private$raw.info); names(has.global) = binary.params

					add.globals = list()
					if (has.column["sample.size"]) {  ## only need global case proportion
						private$configuration$metrics = "sample.size"
						if (has.global["case.proportion"]) add.globals$case.proportion = private$raw.info$case.proportion
						else if (has.global["no.cases"] && has.global["no.controls"]) add.globals$case.proportion = private$raw.info$no.cases / (private$raw.info$no.cases + private$raw.info$no.controls)
					} else {  ## need global sample.size and case.proportion
						if (has.global["sample.size"] && has.global["case.proportion"]) add.globals = private$raw.info[c("sample.size", "case.proportion")]
						else if (has.global["no.cases"] && has.global["no.controls"]) {
							add.globals$sample.size = private$raw.info$no.cases + private$raw.info$no.controls
							add.globals$case.proportion = private$raw.info$no.cases / add.globals$sample.size
						} else if (has.global["sample.size"] && any(has.global[2:3])) {
							add.globals$sample.size = private$raw.info$sample.size
							add.globals$case.proportion = ifelse(has.global["no.cases"], private$raw.info$no.cases, add.globals$sample.size - private$raw.info$no.controls) / add.globals$sample.size
						}
					}

					add.globals = add.globals[!is.na(unlist(add.globals))]
					if (length(add.globals) == 0) private$error("please provide valid combination of input columns and/or global values for %parameter% ", items.and=binary.params)
					invalid.value = !is.finite(unlist(add.globals)) | unlist(add.globals) <= 0 | (names(add.globals) == "case.proportion" & unlist(add.globals) >= 1)
					if (any(invalid.value)) private$error("global %value% provided for %parameter% ", items.and=names(add.globals)[invalid.value], " %is% not valid")

					for (global in names(add.globals)) private$meta.data$sample.metrics[[global]] = add.globals[[global]]
				}

				if ("prevalence" %in% names(private$raw.info)) private$meta.data$sample.metrics$prevalence = private$raw.info$prevalence
			} else {
				binary.params = c("no.cases", "no.controls", "case.proportion", "prevalence")
				if (any(binary.params %in% names(private$raw.info))) private$error("%parameter% ", items.and=binary.params[binary.params %in% names(private$raw.info)], " %is% not valid for continuous phenotype")

				if (is.null(private$raw.info$sample.size)) {
					if (!private$data.interface$main$has.available("sample.size")) private$error("missing both column and global value for required parameter 'sample.size'")
					private$configuration$metrics = "sample.size"
				} else private$meta.data$sample.metrics$sample.size = private$raw.info$sample.size
			}
		},

		get.printer = function() {
			printer = ObjectPrinter$new("PhenotypeInfo (raw input)")$parameter.list(name=self$name(), file=private$input.file, type=self$get.traits())
			if (self$is("multi")) printer$add.parameter("multi-input", paste0(c(private$meta.data$multi$name, private$meta.data$multi$value), collapse="="))
			if (self$is("localized")) {
				units = ifelse(length(private$meta.data$localized$values) > 0, paste0(private$meta.data$localized$names, "=", private$meta.data$localized$values, collapse=", "), paste(private$meta.data$localized$names, collapse=", "))
				printer$add.parameter("localized units", units)
			}
			printer$add.parameter("sample size", private$meta.data$sample.metrics$sample.size)

			if (self$is("binary")) {
				if (!is.null(private$meta.data$sample.metrics$case.proportion)) printer$add.line("case proportion: ", round(100*private$meta.data$sample.metrics$case.proportion, 2), "%", add.space=F)
				printer$parameter.list(prevalence=private$meta.data$sample.metrics$prevalence)
			}
			return(printer)
		}
	),
	public = list(
		initialize = function(...) {
			super$initialize(...)

			input.type = private$raw.info$input.type
			if (length(input.type) != 1 && input.type %in% c("continuous", "binary")) private$error("no valid input.type specified")
			private$traits[[input.type]] = TRUE

			private$validate.files()
			if ("properties" %in% names(private$raw.info)) private$validate.properties()
			if ("parameters" %in% names(private$raw.info)) private$validate.columns()
			if ("unit" %in% names(private$raw.info)) private$validate.localized()

			private$core.columns()
			private$statistic.columns()
			private$metric.columns()
		},

		abbreviate = function() {return(paste0(super$abbreviate(), " (", private$input.file, ")"))},

		name = function(base.only=F) {
			if (!base.only) return(paste0(c(private$pheno.name, private$meta.data$multi$value, private$meta.data$localized$values), collapse="::"))
			else return(private$pheno.name)
		},

		filename = function(chromosomes=NULL) {return(if(!is.null(chromosomes)) gsub("[CHR]", chromosome, private$input.file, fixed=T) else private$input.file)},

		get.global = function(name) {out = private$meta.data$sample.metrics[[name]]; return(if (!is.null(out)) out else NA)},
		get.meta = function(type=c("all", "sample.metrics", "localized", "multi")) {type = match.arg(type); return(if (type == "all") private$meta.data else private$meta.data[[type]])},

		## return list of available parameter names per category; convert to data.frame with corresponding column name if add.columns=T
		get.configuration = function(add.columns=F) {
			if (add.columns) {
				out = data.frame(category=unlist(lapply(names(private$configuration), function(c) {rep(c, length(private$configuration[[c]]))})), stringsAsFactors=F)
				out$parameter = unlist(private$configuration)
				out$column = unlist(private$data.interface$main$get.map(out$parameter, use="first", drop.missing=F))
				return(out)
			} else return(private$configuration)
		},

		get.data = function(chromosomes=NULL) {
			if (!is.null(chromosomes)) {
				if (!self$by.chromosome()) private$error("cannot load by chromosome for whole genome input")
				chromosomes = validate.chromosomes(chromosomes)
				missing = chromosomes$difference(self$get.chromosomes())
				if (!missing$is.empty()) private$error("cannot load data for %chromosome% ", items.and=missing$to.string())
				interfaces = private$data.interface$chromosomes[chromosomes$get()]
			} else if (is.null(chromosomes)) {
				if (self$by.chromosome()) interfaces = private$data.interface$chromosomes[!unlist(lapply(private$data.interface$chromosomes, is.null))]
				interfaces = list(private$data.interface$main)
			}

			if (length(interfaces) == 0) private$error("no data to load")
			for (i in seq_along(interfaces)) {
				curr.input = interfaces[[i]]$load.input(input.type = "summary statistics", error.tag = paste0("phenotype: ", self$name(), "; file: ", interfaces[[i]]$filename()))
				curr.data = curr.input$get.subset(unlist(private$configuration))
				if (i > 1) {
					if (!identical(names(input.data), names(curr.data))) private$error("cannot merge chromosome input data, headers are incompatible", filename=private$input.file)
					else input.data = rbind(input.data, curr.data)
				} else input.data = curr.data
			}

			return(input.data)
		},

		by.chromosome = function() {return(grepl("[CHR]", private$input.file, fixed=T))},
		get.chromosomes = function() {return(ChromosomeList$new(which(sapply(1:23, function(i) {length(private$data.interface$chromosomes) >= i && !is.null(private$data.interface$chromosomes[[i]])}))))},

		## return PhenotypeInfo object with localized units set to specific value; also used for setting multi-input type
		instantiate.units = function(...) {values = lapply(flatten.arglist(...), as.character)
			if (length(values) > 0) {
				if (!any(self$is("localized", "multi"))) private$error("cannot instantiate PhenotypeInfo object with neither 'localized' nor 'multi' type")
				if (any(sapply(values, length) != 1)) private$error("cannot instantiate PhenotypeInfo object, specification is invalid")
				IF.ANY(unknown=!(names(values) %in% c(private$meta.data$multi$name, private$meta.data$localized$names)), THEN=private$error("%unit% ", items.and=names(values)[unknown], " in PhenotypeObject %is% unknown"))

				out = self$clone()
				if (self$is("multi") && private$meta.data$multi$name %in% names(values)) {
					if (!(values[[private$meta.data$multi$name]] %in% names(private$meta.data$multi$columns))) private$error("unit value '", values[[private$meta.data$multi$name]], "' in PhenotypeInfo objecy for multi-input type '", private$meta.data$multi$name, "' is unknown")
					out$.__enclos_env__$private$meta.data$multi$value = values[[private$meta.data$multi$name]]
					values[[private$meta.data$multi$name]] = NULL
				}

				if (length(values) > 0) {
					if (self$is("multi") && is.null(out$.__enclos_env__$private$meta.data$multi$value)) private$error("cannot instantiate PhenotypeInfo object, no value specified for multi-input type '", private$meta.data$multi$name, "'")
					IF.ANY(missing=!(private$meta.data$localized$names %in% names(values)), THEN=private$error("missing unit %value% for localized %unit% ", items.and=private$meta.data$localized$names[missing], " in PhenotypeInfo object"))
					out$.__enclos_env__$private$meta.data$localized$values = values[private$meta.data$localized$names]
				}
				return(out)
			} else return(self)
		}
	)
)


CompositePhenotype = R6::R6Class("CompositePhenotype",
	inherit = PhenotypeInfo,
	private = list(
		input = list(phenotypes = NULL, weights = NULL),

		get.printer = function() {
			label = paste0("PhenotypeInfo (", ifelse(self$is("meta"), "meta", "composite"), ")")
			printer = ObjectPrinter$new(label)$parameter.list(name=self$name(), type=self$get.traits())
			printer$add.list("no. phenotypes", lapply(private$input$phenotypes, function(ph) {ph$abbreviate()}), values.only=T)
			return(printer)
		}
	),
	public = list(
		initialize = function(name, phenotypes, weights, is.meta=F) {
			super$initialize(phenotype=name)
			private$traits[c("composite", if (is.meta) "meta")] = TRUE
			private$input = list(phenotypes = phenotypes, weights = weights)
		},

		get.input = function() {return(private$input)}
	)
)


## instance.id can have multiple components, name will print them separated by ::
SimulatedPhenotype = R6::R6Class("SimulatedPhenotype",
	inherit = PhenotypeInfo,
	private = list(
		instance.id = NULL,  ## eg. iteration ID
		parameters = list(sample.size = NA, heritability = NA),  ## simulation parameters, others may be included dependending on model

		get.printer = function() {
			label = paste0("PhenotypeInfo (simulated)")
			printer = ObjectPrinter$new(label)$parameter.list(name=self$name(), type=self$get.traits())
			printer$add.list("parameters", private$parameters)
			return(printer)
		}
	),
	public = list(
		initialize = function(name, ..., is.binary=F) {
			super$initialize(phenotype=name)
			private$traits[c("simulated", ifelse(is.binary, "binary", "continuous"))] = TRUE
			private$parameters = list.merge(private$parameters, ...)
		},

		name = function(base.only=F) {return(ifelse(!base.only, paste0(c(private$pheno.name, private$instance.id), collapse="::"), private$pheno.name))},

		get.global = function(name) {out = private$parameters[[name]]; return(if (!is.null(out)) out else NA)},
		get.parameters = function() {return(private$parameters)},

		get.id = function() {return(private$instance.id)},
		set.id = function(id) {
			out = self$clone()
			out$.__enclos_env__$private$instance.id = if (length(id) > 0 && !any(is.na(id))) id
			return(out)
		}
	)
)








