

has.type = function(object, types, R6=T, mode=c("all", "any")) {
	if (match.arg(mode) == "any") return(any(types %in% class(object)) && (!R6 || "R6" %in% class(object)))
	else return(all(c(types, if (R6) "R6") %in% class(object)))
}


## input is named arguments (or named list) of format [variable]=[types] (can be multiple types), the variable is checked to see if it has at least one of those types
## to require multiple types, include additional entries for the same variable
## variables are pulled from the environment where check.types() is called; if called (directly) inside a function this will be the internal environment of the function, including the input arguments of that function
## .class can either be a string, or an object with class attribute (using first element); if NULL, will be filled in if a 'self' variable with class R6 is found
## .func is filled in automatically if NULL
check.types = function(..., .func=NULL, .class=NULL, .allow.null=F) {
	env = parent.frame(); args = list(...); label = NULL

	if (is.null(.class) && exists("self", envir=env)) {
		self = get("self", envir=env)
		if (has.type(self, "R6")) .class = class(self)[1]
	} else if (!is.null(.class)) .class = ifelse(is.character(.class), .class[1], class(.class)[1])

	if (is.null(.func)) {
		parent.func = as.list(sys.call(-1))[[1]]
		if (is.call(parent.func)) parent.func = unlist(tail(as.list(parent.func), 1))
		.func = as.character(parent.func)
		if (.func == "initialize") .func = NULL
	}

	if (!is.null(.func)) label = paste0("in function ", if (!is.null(.class)) paste0(.class, "$", collapse=""), .func)
	else if (!is.null(.class)) label = paste0("initializing ", .class, " object")

	if (!is.null(names(args))) {
		for (i in which(names(args) != "")) {
			obj = names(args)[i]
			if (!(obj %in% names(env))) fatal.error("argument '", obj, "' does not exist in input", .label=label)
			if (!has.type(env[[obj]], args[[i]], R6=F, mode="any") && (!.allow.null || !is.null(env[[obj]]))) fatal.error("argument '", obj, "' must have %type% ", items.or=args[[i]], .label=label)
		}
	}
}


## force elements of numeric vector into specific bounds
truncate.values = function(v, min=NULL, max=NULL) {
	if (!is.null(min)) v[!is.na(v) & v < min] = min
	if (!is.null(max)) v[!is.na(v) & v > max] = max
	return(v)
}


## convert input vector to integer values, setting both non-numeric values as well as non-integer numbers to NA
integer.value = function(v) {
	v = suppressWarnings(as.numeric(v))
	v = ifelse(is.infinite(v), v, ifelse(as.integer(v) == v, as.integer(v), NA))
	return(v)
}


## convert vector of integer values to list of disjoint ranges (for sorted values, ignoring non-integer values)
## if to.string=T, render as single string of comma-separated ranges
integer.ranges = function(..., to.string=F) {
	v = sort(unique(integer.value(c(...)))); v = v[!is.na(v)]

	ranges = list()
	if (length(v) > 0) {
		ranges[[1]] = v[1]
		if (length(v) > 1) {
			for (i in 2:length(v)) {
				if (v[i] - v[i-1] == 1) ranges[[length(ranges)]] = c(ranges[[length(ranges)]], v[i])
				else ranges = list.append(ranges, v[i])
			}
		}
		ranges = lapply(ranges, range)
	}

	if (to.string) {
		collapse = function(r) {if (r[1] < 0) r = paste0("(", r, ")") ; return(paste0(r[1], "-", r[2]))}
		ranges = unlist(lapply(ranges, function(r) {if (r[2] - r[1] > 1) collapse(r) else if (r[2] == r[1]) r[1] else r}))
		return(paste(ranges, collapse=", "))
	} else return(ranges)
}

####################


## ... can contain any number of (named) arguments
## arguments evaluate to TRUE if they are boolean values/vectors and contain any TRUE values, or are non-boolean and have a length greater than 0
## if any of the ... arguments evaluate to TRUE, the expression in 'THEN' is executed in the calling environment, attaching any named arguments in ...
## eg. IF.ANY(missing=!(vars %in% var.list), THEN=private$error("missing %variable% ", items.and=var.list[missing]))
## - this will use the private$error function and var.list vector from the calling environment, and the missing vector constructed as function argument
IF.ANY = function(..., THEN) {
	args = list(...)
	execute =	any(unlist(lapply(args, function(arg) {ifelse(is.logical(arg), any(arg), length(arg) > 0)})))
	if (execute) {
		scope = as.environment(clean.list(args)); parent.env(scope) = sys.frame(sys.parent())
		eval(substitute(THEN), envir=scope)
	}
}



## returns named list of input values: unnamed list inputs are recursively flattened to extract elements
## output is a list containing all named elements, and all unnamed non-list elements
## depending on settings, this is subsequently filtered to remove: NA elements, NULL/empty (length == 0) elements, unnamed elements
flatten.arglist = function(..., filter.NA=F, filter.empty=T, filter.unnamed=T) {
	input = list(...); output = list()
	while (length(input) > 0) {
		is.named = !is.null(names(input)[1]) && names(input)[1] != ""
		if (!is.list(input[[1]]) || is.named) {
			output[length(output)+1] = list(input[[1]])
			names(output)[length(output)] = ifelse(is.named, names(input)[1], "")
			input = input[-1]
		} else input = c(input[[1]], input[-1])
	}

	if (length(output) > 0) {
		names(output)[is.na(names(output))] = ""
		drop = rep(F, length(output))
		if (filter.unnamed) drop = drop | names(output) == ""
		if (filter.empty) drop = drop | sapply(output, length) == 0
		if (filter.NA) drop = drop | sapply(output, function(val) {length(val) == 1 && is.na(val)})
		output = output[!drop]
	}
	return(output)
}



##########################

## shorthand to retrieve default settings from globals object
default.settings = function(source, ...) {return(globals$default.settings$get(source, ...))}

## shorthand to retrieve default settings without additional types, and passing ... as override values to default before returning
initialize.settings = function(source, ...) {return(default.settings(source)$set(...))}



## validation helper object for parameter settings stored in objects
## parameters are initialized to default values (named list), and overridden by additional (named) values provided; only parameter names provided in default.values are recognized
## - partial parameter names are matched against default, if it uniquely matched the start of a single parameter
## - multiple values in default or override values are not permitted, except for default of string parameter (treated as enum) or if the 'vector' attribute is explicitly set and is TRUE
## if other (unnamed) ObjectSettings objects are provided, unknown parameters are ignored rather than raising an error
## type of parameter (enum, character, integer, numeric, logical) is automatically derived from default parameter; type can be specified explicitly setting the 'type' attribute on the input value
## - numeric parameters are interpreted as integer if default value is a finite whole number
## - parameter with multiple string values is interpreted as enum, setting first value as default
## minimum and maximum (inclusive) values can be set on numeric/integer parameters via the 'min', 'max' and 'range' attributes
## if default.values is already an ObjectSettings object, returns a modified copy rather than modifying it in place
validate.settings = function(source, default.values, ...) {
	if (has.type(default.values, "ObjectSettings")) return(default.values$set(..., .in.place=F))
	else return(ObjectSettings$new(source, default.values, ...))
}


## STATUS: volatile, shallow
ObjectSettings = R6::R6Class("ObjectSettings",
	private = list(
		source.label = NULL,
		parameters = list(),  ## named list of values for parameters, initialized to default
		param.types = list(),  ## named list of parameter types for parameters: numeric, integer, logical, character, enumeration; can have second element set to "vector", to allow multiple values
		param.range = list(),  ## named list of allowed range of values for parameter, if applicable; list with 'min' and 'max' values for numeric/integer, vector of possible values for enumeration

		error = function(..., .fatal=F) {ifelse(.fatal, fatal.error, input.error)(..., .label=private$source.label)},

		set.defaults = function(default.values) {
			for (i in seq_along(default.values)) {
				param = names(default.values)[i]; value = default.values[[i]]
				if (length(value) > 0 && !any(is.na(value))) {
					if (is.null(attr(value, "type"))) {
						if (is.numeric(value) || is.integer(value)) type = ifelse(is.integer(value) || (all(is.finite(value)) && all(as.integer(value) == value)), "integer", "numeric")
						else if (is.logical(value)) type = "logical"
						else if (is.character(value)) type = ifelse(length(value) > 1, "enumeration", "character")
						else type = class(value)[1]
					} else type = attr(value, "type")[1]

					if (!(type %in% c("numeric", "integer", "logical", "character", "enumeration"))) private$error("parameter '", param, "' has unsupported type '", type, "'", .fatal=T)

					if (!is.null(attr(value, "vector")) && attr(value, "vector")) {
						if (type == "enumeration") private$error("parameter '", param, "' of type 'enumeration' cannot have 'vector' attribute", .fatal=T)
						private$param.types[[param]] = c(type, "vector")
					} else private$param.types[[param]] = type

					if (type != "enumeration") {
						private$parameters[[param]] = value
						if (type == "numeric" || type == "integer") {
							if (!is.null(attr(value, "range"))) range = named.list(attr(value, "range")[1:2], names=c("min", "max"))
							else range = list(min=ifelse(!is.null(attr(value, "min")), attr(value, "min"), -Inf), max=ifelse(!is.null(attr(value, "max")), attr(value, "max"), Inf))
							private$param.range[[param]] = range
						}
					} else {
						private$parameters[[param]] = value[1]
						private$param.range[[param]] = value
					}
				} else private$error("parameter '", param, "' has missing or NULL default value", .fatal=T)
			}
			if (length(private$parameters) > 0) {
				multiple = sapply(private$parameters, length) > 1 & !sapply(private$param.types, function(types) {"vector" %in% types})
				if (any(multiple)) private$error("multiple default values provided for non-vector %parameter% ", items.and=unique(names(private$parameters)[multiple]))
			}
		},

		set.value = function(param, value) {
			if (length(value) > 1 && !("vector" %in% private$param.types[[param]])) private$error("value for non-vector parameter '", param, "' has length greater than one")
			if (is.null(value) || any(is.na(value))) private$error("value for parameter '", param, "' is NULL or contains missing values")

			type = private$param.types[[param]][1]
			if (type == "character" || type == "enumeration") {
				value = try(as.character(value), silent=T)
				if (class(value) == "try-error") private$error("value for parameter '", param, "' cannot be converted to a string")
				if (type == "enumeration") {
					value = try(match.arg(value, private$param.range[[param]]), silent=T)
					if (class(value) == "try-error") private$error("value for parameter '", param, "' should be one of ", items.or=private$param.range[[param]])
				}
			} else if (type == "numeric" || type == "integer") {
				value = suppressWarnings(as.numeric(value))
				if (type == "integer" && !any(is.na(value))) value = ifelse(is.infinite(value), value, ifelse(as.integer(value) == value, as.integer(value), NA))
				if (any(is.na(value))) private$error("value for parameter '", param, "' is not ", ifelse(type == "numeric", "a number", "an integer"))
				if (any(value < private$param.range[[param]]$min)) private$error("value for parameter '", param, "' should be equal to or larger than ", private$param.range[[param]]$min)
				if (any(value > private$param.range[[param]]$max)) private$error("value for parameter '", param, "' should be equal to or smaller than ", private$param.range[[param]]$max)
			} else if (type == "logical") {
				if (!is.logical(value)) private$error("value for parameter '", param, "' should be TRUE or FALSE")
			} else fatal.error("parameter type '", type, "' not implemented")

			private$parameters[[param]] = value
		}
	),
	public = list(
		initialize = function(source, default.values, ...) {
			private$source.label = ifelse(!is.character(source), paste(class(source)[1], "object"), source)
			private$set.defaults(default.values)
			self$set(...)
		},

		print = function(...) {
			printer = ObjectPrinter$new("ObjectSettings")$parameter.list(source=private$source.label)
			params = private$parameters
			if (length(params) > 0) names(params) = paste0(names(params), " (", sapply(private$param.types[names(params)], paste, collapse=", "), ")")
			printer$add.list("no. parameters", params)
			cat(printer$to.string())
		},

		get = function(parameter) {
			if (!(parameter %in% names(private$parameters))) private$error("trying to retrieve unknown parameter '", parameter, "'")
			return(private$parameters[[parameter]])
		},
		get.all = function() {return(private$parameters)},

		set.source = function(source, .in.place=T) {
			if (.in.place) private$source.label = ifelse(!is.character(source), paste(class(source)[1], "object"), source)
			else self$clone()$set.source(source, .in.place=T)
		},

		## input (combination of) parameter=value named arguments and named lists with parameter=value elements
		## if .in.place=F, will return a copy of itself with modified values, rather than modifying in place
		set = function(..., .in.place=T) {input = list(...)
			if (.in.place) {
				is.settings = as.logical(unlist(lapply(input, function(i) {has.type(i, "ObjectSettings")})))
				for (i in which(is.settings)) {
					add = input[[i]]$get.all()
					input[[i]] = add[names(add) %in% names(private$parameters)]
				}
				input = flatten.arglist(input, filter.empty=F)

				if (length(input) > 0) {
					names(input) = complete.names(names(input), names(private$parameters))
					IF.ANY(unknown = !(names(input) %in% names(private$parameters)), THEN=private$error("unknown %parameter% ", items.and=unique(names(input)[unknown])))
					IF.ANY(multiple = sapply(input, length) > 1 & !sapply(private$param.types[names(input)], function(types) {"vector" %in% types}), THEN=private$error("multiple values provided for non-vector %parameter% ", items.and=unique(names(input)[multiple])))
					IF.ANY(empty = sapply(input, length) == 0, THEN=private$error("%value% for %parameter% ", items.and=unique(names(input)[empty]), " %is% NULL or empty vector"))
					IF.ANY(missing = sapply(input, function(val) {any(is.na(val))}), THEN=private$error("%value% for %parameter% ", items.and=unique(names(input)[missing]), " contains missing values"))

					for (i in 1:length(input)) private$set.value(names(input)[i], input[[i]])   ## if duplicate names, last instance wins
				}
				invisible(self)
			} else return(self$clone()$set(..., .in.place=T))
		},

		## merge in additional settings objects; unlike set(), will add new parameters not already in ObjectSettings object
		merge = function(..., .in.place=T) {args = flatten.arglist(..., filter.unnamed=F)
			if (length(args) > 0) {
				if (!all(sapply(args, has.type, "ObjectSettings"))) private$error("cannot merge settings, all input must be of type ObjectSettings")
				if (.in.place) {
					for (other in args) {
						other = other$.__enclos_env__$private
						for (var in c("parameters", "param.types", "param.range") ) private[[var]] = list.merge(private[[var]], other[[var]])
					}
					invisible(self)
				} else return(self$clone()$merge(..., .in.place=T))
			} else invisible(self)
		},

		## run parameter value through settings object, return existing value of input value is NULL and null.default=T
		validate = function(parameter, value, null.default=T) {
			if (!is.null(value) || !null.default) return(self$set(named.list(parameter, value), .in.place=F)$get(parameter))
			else return(self$get(parameter))
		}
	)
)


## storage object for default settings, by class name
## for add(), character inputs with IMPORT attribute set to TRUE denote class names automatically included when calling get() for that class (inserted before class defaults)
## get() will retrieve and merge all available settings objects matching class of specified source, plus additional specified in ... (class names, or class of additional objects)
DefaultSettings = R6::R6Class("DefaultSettings",
	private = list(default.settings = list(), auto.import = list()),
	public = list(
		add = function(class.name, ...) {
			import = unlist(lapply(list(...), function(arg) {if (is.character(arg) && !is.null(attr(arg, "IMPORT")) && attr(arg, "IMPORT")) arg}))
			if (length(import) > 0) private$auto.import[[class.name]] = import
			private$default.settings[[class.name]] = validate.settings(class.name, flatten.arglist(...))
		},

		get = function(source, ...) {
			settings = ObjectSettings$new(source, list())

			process.list = unique(unlist(lapply(list(source, ...), function(t) {if (!is.character(t)) rev(class(t)) else t})))
			process.list = rev(process.list[process.list %in% names(private$default.settings)])

			include = c()
			while (length(process.list) > 0) {
				include = c(process.list[1], include); import = private$auto.import[[process.list[1]]]
				import = import[import %in% names(private$default.settings) & !(import %in% include)]
				process.list = c(rev(import), process.list[-1])
			}

			if (length(include) > 0) settings$merge(private$default.settings[include])
			return(settings)
		}
	)
)



################################



## completes elements in 'targets' to match element in 'range' (ie. exact match, or matches start of single element in 'range')
## if no match or multiple matches, returns original value, or NA if failed.NA=T
complete.names = function(targets, range, failed.NA=F) {
	exact = targets %in% range
	for (i in which(!exact)) {
		partial = targets[i] == substr(range, 1, nchar(targets[i]))
		if (sum(partial) == 1) targets[i] = range[partial]
		else if (failed.NA) targets[i] = NA
	}
	return(targets)
}


## performs partial match of targets to range, returns list with valid matches, as well as ambiguous and unknown targets
## 'matches' output value is in same order as 'targets' input, and same length if none failed and all entries in 'targets' are unique and no wildcards were expanded
## wildcard * matches any number of characters, wildcard ? matches a single character
partial.match = function(targets, range, expand.wildcard=F) {
	matches = list(); is.wild = c()
	for (t in targets) {
		if (!any(range == t)) {
			curr = unique(grep(t, range, value=T, fixed=T))
			matches[[t]] = curr[substr(curr, 1, nchar(t)) == t]
		} else matches[[t]] = t

		if (expand.wildcard && length(matches[[t]]) == 0 && grepl("\\*|\\?", t)) {
			pattern = gsub("\\?", ".", gsub("\\*", ".*", gsub("\\.", "\\\\.", t))) ## escape periods, convert wildcards to regular expression
			pattern = paste0("^", pattern, "$")
			res = try(grep(pattern, range, value=T), silent=T)
			if (class(res) != "try-error")	matches[[t]] = res
			if (length(matches[[t]]) > 0) is.wild = c(is.wild, t)
		}
	}

	count = sapply(matches, length)
	if (length(is.wild) > 0) count[count > 1 & names(matches) %in% is.wild] = 1
	return(list(
		matches = as.vector(unlist(matches[count == 1])),
		ambiguous = targets[count > 1],
		unknown = targets[count == 0],
		failure = any(count != 1)
	))
}


## scrub whitespace around commas and equality signs, and split by whitespace / semi-colons, removing empty elements
preprocess.specifier = function(input.str) {
	input.str = gsub(",+", ",", gsub("\\s*,\\s*", ",", gsub("\\s*=\\s*", "=", trimws(input.str))))   ## remove whitespace around commas and equality signs, and collapse multiple commas
	elems = strsplit(trimws(input.str), "(\\s*;\\s*)|(\\s+)")[[1]]   ## split by whitespace / semicolons
	elems = gsub(",*=,*", "=", gsub("(^,)|(,$)", "", elems))   ## remove extraneous commas
	return(elems[elems != ""])
}

## parse input.string of the format "name1=value1,value2,... name2=value3,value4,...; name3=value5 value6" etc.
## whitespace within primary name/value elements around commas and equality signs is ignored
## primary are separated by whitespace and/or semi-colons
## if infer.names=T, unnamed elements like 'value6' are treated as 'value6=value6'
parse.value.string = function(input.str, allow.unnamed=T, infer.names=allow.unnamed) {
	raw = preprocess.specifier(input.str)
	pairs = strsplit(raw, "="); counts = sapply(pairs, length)

	invalid = counts > 2 | grepl("(^=)|(=$)", raw)  ## more than one '=', or '=' at beginning or end of element
	if (any(invalid)) input.error("invalid specifier %token% ", items.and=raw[invalid])

	output = strsplit(sapply(pairs, tail, 1), ",")
	names(output) = sapply(pairs, function(spec) {ifelse(length(spec) == 2, spec[1], "")})

	if (infer.names) {
		infer = names(output) == "" & sapply(output, length) == 1
		if (any(infer)) names(output)[infer] = unlist(output[infer])
	}

	if (!allow.unnamed && any(names(output) == "")) input.error("specification string '", input.str, "' contains unnamed ", if (infer.names) "multi-value ", "elements")
	return(output)
}



## check list of (partial) phenotype names against list of available phenotypes
## expands * wildcard to match multiple characters if not already matched as literal *, and same for expanding ? wildcard to match single character
## throws error if any are ambiguous or unmatched, returns vector with corresponding full phenotype names otherwise (removing any duplicates)
## returns NULL if phenotypes argument is NULL and character(0) if phenotypes argument is empty
check.phenotypes = function(phenotypes, available, label=NULL, discard.ambiguous=F, discard.unknown=F) {
	if (length(phenotypes) > 0) {
		index = partial.match(phenotypes, available, expand.wildcard=T)

		if (!is.null(label)) label = paste0(" in ", label)
		if (!discard.ambiguous && length(index$ambiguous) > 0) input.error("ambiguous partial phenotype %name% ", items=index$ambiguous, " %match[-es]% multiple phenotypes", label)
		if (!discard.unknown && length(index$unknown) > 0) input.error("%phenotype% ", items.and=index$unknown, " %is% not present", label)
		return(unique(index$matches))
	} else {
		if (is.null(phenotypes)) return(NULL)
		else return(character(0))
	}
}


#########################


## if resolve.chr=T and filename contains [CHR] token, considered to exist if at least one chromosome-specific file matching the filename template is found
check.files.exist = function(..., resolve.chr=F) {
	file.names = c(...)
	if (length(file.names) > 0) {
		exist = !is.na(file.names) & file.exists(file.names)
		if (resolve.chr) {
			tpl = !is.na(file.names) & grepl("[CHR]", file.names, fixed=T)
			if (any(tpl))	exist[tpl] = unlist(lapply(file.names[tpl], function(f) {nrow(chromosome.files(f)) > 0}))
		}
		if (!all(exist)) input.error("missing input %file% ", items.and=file.names[!exist], "; please ensure correct file names and paths have been provided")
	}
	invisible(TRUE)
}

check.dirs.exist = function(...) {
	dirs = c(...)
	if (length(dirs) > 0) {
		valid = !is.na(dirs) & dir.exists(dirs)
		if (!all(valid)) input.error("missing or invalid %directory% ", items=dirs[!valid], "; please ensure correct paths have been provided")
	}
	invisible(TRUE)
}



#########################



chromosome.files = function(file.tpl, chromosomes="all", prune.missing=T) {
	if (grepl("[CHR]", file.tpl, fixed=T)) {
		files = sapply(c(1:23, "X", "x"), function(c) {gsub("[CHR]", c, file.tpl, fixed=T)})
		files[!file.exists(files)] = NA
		chr.x = files[23:25]
		if (any(!is.na(chr.x))) files[23] = chr.x[!is.na(chr.x)][1]
		out = data.frame(chromosome=1:23, file=files[1:23], stringsAsFactors=F)

		chromosomes = validate.chromosomes(chromosomes)
		if (!chromosomes$has.all()) out = out[out$chromosome %in% chromosomes$get(),]
		if (prune.missing) out = out[!is.na(out$file),]
		return(out)
	} else return(NULL)
}

## if condense=T, validate input and convert to ChromosomeList object
## otherwise, return vector of same length as 'chromosomes'; if any chromosomes are invalid, if mode="filter" these are set to NA, otherwise an error/warning is generated
validate.chromosomes = function(chromosomes, condense=T, mode=c("error", "warning", "filter")) {
	if (condense) return(if (!has.type(chromosomes, "ChromosomeList")) ChromosomeList$new(chromosomes) else chromosomes)
	else return(ChromosomeList$new("all")$validate(chromosomes, mode=match.arg(mode)))
}


## wrapper object for list of chromosome codes; allows for easy validation, printing, and comparison/combination of chromosome lists
## uses numeric coding (X is converted to 23); initialize with chromosomes="none" for empty list
## member functions with 'chromosomes' argument can take either a vector of chromosome codes, or another ChromosomeList object
ChromosomeList = R6::R6Class("ChromosomeList",
	private = list(
		chromosomes = c(),

		set.chromosomes = function(chromosomes) {
			if (!is.null(chromosomes) && !any(chromosomes == "all")) {
				chromosomes[tolower(chromosomes) == "x"] = 23
				unknown = !(chromosomes %in% 1:23)
				if (any(unknown)) input.error("unknown chromosome %code% ", items.and=unique(chromosomes[unknown]))
				private$chromosomes = sort(unique(as.numeric(chromosomes)))
			} else private$chromosomes = 1:23
		},

		to.list = function(chromosomes) {return(if (has.type(chromosomes, "ChromosomeList")) chromosomes else ChromosomeList$new(chromosomes))}
	),
	public = list(
		initialize = function(chromosomes="all") {if (length(chromosomes) > 0 && chromosomes[1] != "none") private$set.chromosomes(chromosomes)},
		print = function(...) {cat(self$to.string())},

		to.string = function() {
			if (!self$has.all()) {
				if (self$size() > 0) {
					out = if (self$size(auto.only=T) > 0) integer.ranges(private$chromosomes[private$chromosomes != 23], to.string=T)
					if (self$has.x()) out = paste(c(out, "X"), collapse=", ")
					return(add.plural(out, is.plural = (self$size() > 1)))
				} else return("[none]")
			} else return(add.plural("1-22, X", is.plural=TRUE))
		},

		get = function() {return(private$chromosomes)},
		size = function(auto.only=F) {out = length(private$chromosomes); return(ifelse(auto.only, out-self$has.x(), out))},
		is.empty = function() {return(self$size() == 0)},

		has.x = function() {return(23 %in% private$chromosomes)},
		has.all = function() {return(length(private$chromosomes) == 23 && all(private$chromosomes == 1:23))},

		equals = function(chromosomes) {chromosomes = private$to.list(chromosomes); return(self$size() == chromosomes$size() && all(private$chromosomes == chromosomes$get()))},

		contains = function(chromosomes, aggregate=T) {
			if (is.character(chromosomes) && length(chromosomes) == 1 && chromosomes == "all") return(self$has.all())
			if (has.type(chromosomes, "ChromosomeList")) chromosomes = chromosomes$get()
			return(if (aggregate) all(chromosomes %in% private$chromosomes) else chromosomes %in% private$chromosomes)
		},

		## difference() returns chromosomes in current ChromosomeList not contained in other
		difference = function(chromosomes) {chromosomes = private$to.list(chromosomes);	return(ChromosomeList$new(setdiff(private$chromosomes, chromosomes$get())))},
		union = function(chromosomes) {chromosomes = private$to.list(chromosomes);	return(ChromosomeList$new(union(private$chromosomes, chromosomes$get())))},
		intersection = function(chromosomes) {chromosomes = private$to.list(chromosomes);	return(ChromosomeList$new(intersect(private$chromosomes, chromosomes$get())))},

		## validate chromosome codes and check against current set
		validate = function(chr.vec, mode=c("filter", "warning", "error")) {mode = match.arg(mode)
			chr.vec[chr.vec == "x" | chr.vec == "X"] = 23
			invalid = !is.na(chr.vec) & !(chr.vec %in% private$chromosomes)
			if (mode != "filter" && any(invalid)) ifelse(mode == "warning", throw.warning, input.error)("unknown or unavailable chromosome %code% ", items.and=unique(chr.vec[invalid]), " in input")
			chr.vec[invalid] = NA
			return(as.integer(chr.vec))
		}
	)
)








