## InputInterface provides an interface that maps a set of variables (eg. column names) onto internal parameters
## the aim is to map potentially heterogeneous data input onto a fixed set of parameters, resolving any conflicts and duplications

## the FileInputInterface class (and extensions) below provides an InputInterface specifically for input files, mapping parameters to columns in a file
## the InputData class (and extensions) below provide mapped access to a data source (data.frame or input file) using an InputInterface
## - this allows data columns to be accessed directly via their parameter names


#########################


## the 'input.names' argument specifies the available variable names
## the 'parameter.index' argument defines the mapping of parameters onto input names (see eg. the *.header.index lists in globals.R)
## - the parameter index should be a named list: names are the parameter names, values are vectors of possible input names that map to that parameter
## - mapping is case insensitive, and every input name should occur only once, for a single parameter
## - parameter names are included as possible input names, unless include.param.names=F

## the get.available(), get.missing() and get.duplicates() functions can be used to determine whether specific parameters are mapped to any input names, and whether any are mapped to multiple input names
## the get.map() function returns a map of parameters onto available input names, by default selecting the first mapped input name for each parameter and trimming away unmapped parameters
InputInterface = R6::R6Class("InputInterface",
	private = list(
		input.names = NULL,  ## ordered vector of variable names to map
		input.type = NULL,  ## label for the input units (eg. 'variable', 'column')
		parameter.index = list(),  ## named (parameters) list: vectors of possible input names per parameter
		parameter.map = list(),  ## named (parameters) list: vectors of actual input names available for each parameter (stored in original case)
		match.counts = NULL,  ## named (parameters) list: counts of number of mapped input names per parameter

		create.map = function() {
			IF.ANY(duplicates = duplicated(unlist(private$parameter.index)), THEN = input.error("invalid index, %", private$input.type, "% ", items.and=unique(unlist(private$parameter.index)[duplicates]), " %is% %[+each]% mapped to multiple parameters"))

			private$parameter.map = list()
			for (h in names(private$parameter.index)) private$parameter.map[[h]] = unique(private$input.names[tolower(private$input.names) %in% private$parameter.index[[h]]])

			private$match.counts = unlist(lapply(private$parameter.map, length))
		}
	),
	public = list(
		initialize = function(input.names, parameter.index, input.type="variable", include.param.names=T) {
			private$input.names = input.names; private$input.type = input.type

			if (include.param.names) {for (h in names(parameter.index)) parameter.index[[h]] = c(parameter.index[[h]], h)}
			private$parameter.index = clean.list(lapply(parameter.index, function(cols) {unique(tolower(cols))}))

			private$create.map()
		},

		get.variables = function() {return(private$input.names)},
		get.parameters = function() {return(names(private$parameter.index))},

		get.available = function() {return(self$get.parameters()[private$match.counts > 0])},
		get.missing = function() {return(self$get.parameters()[private$match.counts == 0])},
		get.duplicates = function() {return(self$get.parameters()[private$match.counts > 1])},

		has.parameters = function(parameters, aggregate=F) {status = parameters %in% self$get.parameters(); return(if (aggregate) all(status) else status)},
		has.available = function(parameters, aggregate=F) {status = parameters %in% self$get.available(); return(if (aggregate) all(status) else status)},

		## return map to first/last/all mapped input elements, trimming away unmapped parameters unless trim.empty=F
		## if parameters=NULL, return full map of all parameters, otherwise subset to selected parameters
		## if drop.missing=T, unknown parameters are discarded, otherwise they are included with an NA value
		get.map = function(parameters=NULL, use=c("first", "last", "all"), drop.missing=T, trim.empty=T) {use = match.arg(use)
			if (use == "all") map = private$parameter.map
			else map = lapply(private$parameter.map, ifelse(use == "first", head, tail), 1)

			if (trim.empty) map = map[private$match.counts > 0]
			if (!is.null(parameters)) {
				IF.ANY(unknown = !(parameters %in% self$get.parameters()), THEN = input.error("unknown %parameter% ", items.and=parameters[unknown], " in get.map() function of InputInterface"))
				map = map[names(map) %in% parameters]
				if (!drop.missing) map[parameters[!(parameters %in% names(map))]] = NA
				map = map[order(match(names(map), parameters))]
			}
			return(map)
		},

		## return data.frame with partial matched for each of the specified parameters, with columns parameter, input.name, label
		## variables are matched if their name is of the form [INPUT_NAME]_[LABEL] (if sep="_")
		## unless filter="none", column names that are already either in the parameter.index or in the parameter.map are excluded
		map.partial = function(parameters=NULL, sep="_", filter=c("index", "map", "none")) {filter = match.arg(filter)
			if (!is.null(parameters)) IF.ANY(unknown = !(parameters %in% self$get.parameters()), THEN = input.error("unknown %parameter% ", items.and=parameters[unknown], " in partial.map() function of InputInterface"))
			else parameters = self$get.parameters()

			matched = data.frame(parameter=character(0), input.name=character(0), label=character(0))
			for (param in parameters) {
				for (prefix in paste0(private$parameter.index[[param]], sep)) {
					found = grep(paste0("^", prefix, ".+"), private$input.names, ignore.case=T, value=T)
					if (length(found) > 0) matched = rbind(matched, data.frame(parameter=param, input.name=found, label=tolower(gsub(paste0("^", prefix), "", found, ignore.case=T))))
				}
			}
			matched = matched[!duplicated(paste(matched$parameter, matched$label)),]
			if (filter != "none") matched = matched[!(tolower(matched$input.name) %in% unlist(if (filter == "index") private$paramater.index else private$parameter.map)),]
			return(matched)
		},

		## add additional parameter mappings (same format as 'parameter.index'), adding new parameters
		## if append.existing=T, mapping of existing parameters is extended with additional input names; otherwise, they are ignored
		## returns a list of parameters included in the input mapping that already existed
		update.parameters = function(add.index, append.existing=T) {
			add.index = clean.list(add.index); existing = c()
			for (param in names(add.index)) {
				if (param %in% names(private$parameter.index)) {
					existing = c(existing, param)
					if (!append.existing) next
				}
				private$parameter.index[[param]] = unique(c(private$parameter.index[[param]], tolower(add.index[[param]])))
			}
			private$create.map()
			return(existing)
		}
	)
)


## InputInterface for a file header, mapping parameters onto column names; reads in and parses the file header on initialization
## also reads in global parameters in the file above the file header (if any)
FileInputInterface = R6::R6Class("FileInputInterface",
	inherit = InputInterface,
	private = list(
		input.file = NULL,
		file.header = NULL, ## vector of column names
		global.parameters = list(),  ## named list of global parameters read from meta information at top of file

		read.header = function() {undefined.error("read.header", class(self)[1])}
	),
	public = list(
		initialize = function(filename, parameter.index, include.param.names=T) {
			if (check.files.exist(filename)) private$input.file = filename
			private$read.header()
			super$initialize(private$file.header, parameter.index, input.type="column",	include.param.names=include.param.names)
		},

		filename = function() {return(private$input.file)},
		get.globals = function() {return(private$global.parameters)}
	)
)


## specific implementation for fixed-width, text files
## the header line is identified as the first line not starting with a # character (ignoring leading whitespace)
## global parameters above the header are identified as any line of the form "# [NAME]: [VALUES]", ignoring any leading/trailing whitespace or whitespace around the colon
## the [VALUES] component is cleaned up as follows (see also preprocess.specifier function):
## - leading/trailing whitespace is removed, as is any whitespace around equality signs and commas, collapsing consecutive commas into one comma
## - it is then split into individual elements on semi-colons and whitespace
## - for each element, commas adjacent to an equality sign or at the start or end of the element are removed; empty elements are then discarded
## "# settings: values=100, 200, 300; sort = numeric DEBUG" would yield a global parameter 'settings' with value c("values=100,200,300", "sort=numeric", "DEBUG")
TextFileInputInterface = R6::R6Class("TextFileInputInterface",
	inherit = FileInputInterface,
	private = list(
		comment.lines = 0,  ## number of comment lines preceding the header

		read.header = function() {
			comments = NULL; curr.step = 1; increment = 100
			while (length(private$file.header) == 0) {
				current = scan(private$input.file, sep="\n", what="", quiet=T, strip.white=T, blank.lines.skip=F, skip=(curr.step-1)*increment, nlines=increment)
				if (length(current) == 0) break

				input.lines = grep("^[^#]", current) ## match a non-empty line not starting with a # (NB: leading/trailing whitespaces are stripped by scan, so empty lines are "")
				if (length(input.lines) > 0) {
					if (input.lines[1] > 1) {
						comments = c(comments, current[1:(input.lines[1]-1)])
						private$comment.lines = length(comments)
					}
					private$file.header = strsplit(current[input.lines[1]], "\\s+")[[1]]
					break
				} else {
					comments = c(comments, input.lines)
					curr.step = curr.step + 1
				}
			}
			if (length(private$file.header) == 0) input.error("input file ", private$input.file, " contains no header")

			comments = trimws(gsub("^#", "", comments))
			parameters = grep("^\\S+\\s*:\\s*\\S+", comments[comments != ""], value=T)   ## select parameter comments (non-whitespace string, (optional whitespace), colon, (optional whitespace), any amount of non-whitespace

			if (length(parameters) > 0) {
				param.parts = lapply(strsplit(parameters, ":"), function(v) {trimws(c(v[1], paste(v[-1], collapse=":")))})  ## split on first colon
				for (i in 1:length(param.parts)) {
					values = preprocess.specifier(param.parts[[i]][2])
					if (length(values) > 0) private$global.parameters[[param.parts[[i]][1]]] = values
				}
				private$global.parameters = clean.list(private$global.parameters)
			}
		}
	),
	public = list(
		## returns a full TextFileInput object for the input file
		load.input = function(input.type=NULL, error.tag=NULL, defer.loading=F) {return(TextFileInput$new(self, input.type=input.type, error.tag=error.tag, defer.loading=defer.loading))},

		## returns a data.frame with the data in the input file
		load.data = function() {
			data = data.table::fread(private$input.file, data.table=F, showProgress=F, skip=private$comment.lines)
			if (length(private$file.header) != ncol(data)) input.error("invalid header for file ", private$input.file)
			return(data)
		}
	)
)


####################################


## wraps around a data source, and provides access to underlying data via parameter mapping as defined by InputInterface
## base class is not for direct construction, use List/Dataframe/Textfile derived classes
InputData = R6::R6Class("InputData",
	private = list(
		data = NULL, ## data.frame or list
		param.interface = NULL, ## InputInterface object mapping parameter names to data column names
		use.mode = NULL, ## which value to use if multiple, can be 'first' or 'last'
		error.tag = NULL,

		error = function(..., fatal=F) {ifelse(fatal, fatal.error, input.error)(..., .label=private$error.tag)},

		match.parameters = function(..., target, as.bool=F, check.unknown=T) {parameters = c(...)
			unknown = !private$param.interface$has.parameters(parameters)
			if (check.unknown && any(unknown)) private$error("%parameter% ", items.and=parameters[unknown], " %is% unknown")
			if (as.bool)	return(parameters %in% target)
			else return(parameters[parameters %in% target])
		},

		check.required = function(...) {
			missing = private$match.parameters(..., target=private$param.interface$get.missing())
			if (length(missing) > 0) private$error("no input %column% %was% found for required %parameter% ", items.and=missing)
		},

		check.duplicates = function(...) {
			duplicates = private$match.parameters(..., target=private$param.interface$get.duplicates())
			if (length(duplicates) > 0) throw.warning("matched multiple columns for %parameter% ", items.and=duplicates, ", using first matched column %[+for each]%", .label=private$error.tag)
		},

		## return data columns corresponding to parameters specified in ... argument
		## columns are renamed to parameter names, and the first occurrence is used if multiple columns match the same parameter
		extract.parameters = function(...) {
			private$check.required(...)
			private$check.duplicates(...)

			index = private$param.interface$get.map(c(...), use=private$use.mode, drop.missing=T)
			return(add.names(private$data[unlist(index)], names(index)))
		}
	),
	public = list(
		initialize = function(data, interface, use.mode=c("first", "last"), error.tag=NULL) {
			private$data = data; private$param.interface = interface; private$use.mode = match.arg(use.mode)
			if (!is.null(error.tag)) private$error.tag = error.tag;
			IF.ANY(duplicate = duplicated(tolower(interface$get.variables())), THEN = private$error("duplicate column %name% ", items.and=interface$get.variables()[duplicate]))
		},

		get.variables = function() {return(private$param.interface$get.variables())},
		get.parameters = function(mode=c("available", "all")) {return(if (match.arg(mode) == "available") private$param.interface$get.available() else private$param.interface$get.parameters())},

		has.parameters = function(..., aggregate=T) {
			available = private$match.parameters(..., target=private$param.interface$get.available(), as.bool=T, check.unknown=F)
			return(if (aggregate) all(available) else available)
		},

		## get data.frame with data for all available parameters
		get.data = function() {return(private$extract.parameters(private$param.interface$get.available()))},

		## get data.frame with data for specified parameters; parameters named in ... are required, and an error is thrown if any are not available
		get.subset = function(..., optional=NULL) {
			parameters = c(...)
			if (!is.null(optional)) parameters = c(parameters, optional[private$param.interface$has.available(optional)])
			return(private$extract.parameters(parameters))
		},

		## get specific parameter, and return as vector
		get.parameter = function(param.name) {return(private$extract.parameters(param.name)[[1]])}
	)
)


## wraps around a named list (unnamed entries are disregarded)
## entries with NULL values are discarded (default) or converted to NA
ListInput = R6::R6Class("ListInput",
	inherit = InputData,
	public = list(
		initialize = function(data, header.index, drop.null=T, use.mode=c("first", "last"), input.type=NULL, error.tag=NULL) {
			if (is.null(error.tag)) private$error.tag = paste0("in ", if (!is.null(input.type)) paste0(input.type, " "), "list input")
			if (!is.list(data)) private$error("input to ListInput object is not a list", fatal=T)
			data = clean.list(as.list(data), drop.null=drop.null)
			if (!drop.null) data[unlist(lapply(data, is.null))] = NA

			interface = InputInterface$new(names(data), header.index, include.param.names=T)
			super$initialize(data, interface, use.mode=match.arg(use.mode))
		}
	)
)




## generic tabular data base class
TabularData = R6::R6Class("TabularInput", inherit = InputData)


## wraps around input data.frame
DataframeInput = R6::R6Class("DataFrameInput",
	inherit = TabularData,
	public = list(
		initialize = function(data, header.index, use.mode=c("first", "last"), input.type=NULL, error.tag=NULL) {
			if (is.null(error.tag)) private$error.tag = paste0("in ", if (!is.null(input.type)) paste0(input.type, " "), "data.frame input")
			if (!is.data.frame(data)) private$error("input to DataframeInput object is not a data.frame", fatal=T)

			interface = InputInterface$new(names(data), header.index, input.type="column", include.param.names=T)
			super$initialize(data, interface, use.mode=match.arg(use.mode))
		}
	)
)


## wraps around text file; data loading can be deferred until first data request
TextFileInput = R6::R6Class("TextFileInput",
	inherit = TabularData,
	private = list(
		input.file = NULL,

		extract.parameters = function(...) {
			if (is.null(private$data)) self$load.data()
			return(super$extract.parameters(...))
		}
	),
	public = list(
		initialize = function(input, header.index, use.mode=c("first", "last"), input.type=NULL, error.tag=NULL, defer.loading=F) {
			if (has.type(input, "TextFileInputInterface")) interface = input
			else interface = TextFileInputInterface$new(input, header.index, include.param.names=T)
			private$input.file = interface$filename()

			if (is.null(error.tag)) error.tag = paste0("in ", if (!is.null(input.type)) paste0(input.type, " "), "file ", private$input.file)
			super$initialize(NULL, interface, use.mode=match.arg(use.mode), error.tag=error.tag)
			if (!defer.loading) self$load.data()
		},

		get.globals = function() {return(private$param.interface$get.globals())},
		load.data = function() {if (is.null(private$data)) private$data = private$param.interface$load.data()}
	)
)







