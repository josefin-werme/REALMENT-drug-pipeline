## log, warning and error messages all use the same underlying ConditionMessage structure
## - message is concatenated together from ... input, with named arguments starting with 'items' pasted into a comma-separated list of single quoted items (items.and/items.or will add and/or between last two elements)
##   - tokens surrounded by %'s in the concatenated input are pluralized if any 'items' arguments have length > 1 (NB: irregular plural forms must be added to globals$plural.form in globals.R)
##   - existing ConditionMessage and base R error objects (eg. caught errors from tryCatch) can be passed into these functions as well (eg. a caught error can be directly converted to a warning by passing it to throw.warning)
## - special arguments (general):
##   - .prefix: appended in front of the main message (NB: any spacing between prefix and main message must be specified explicitly in .prefix)
##   - .labels: vector or (named) list of values, appended after main message in parentheses as semi-colon separated list; prints as "name: value" if named list item
##   - .lines: vector of lines to print on separate lines below the main message
## - special arguments (only used for errors)
##   - .data: arbitrary list embedded in the ConditionMessage object to transfer additional information
##   - .sub.type: error subtypes, can be used as tryCatch handler name to catch specific errors

## the form of the intendation is specified by the log.indent.string in the globals list in global.R
log.message = function(..., .indent=1) {
	msg = create.message(...)$set(.prefix=paste0(rep(globals$log.indent.string, times=.indent), collapse=""))$get.message()
	if (globals$warning.message.state) msg = paste0("\n", msg)
	message(msg); globals$warning.message.state = F
}

## .tag argument is printed in uppercase in front of warning message; defaults to WARNING or converted error type if not explicitly set
throw.warning = function(..., .tag=NULL) {
	msg = create.message(...)
	if (is.null(.tag) && is.error(msg)) {
		error.type = grep(".*\\.error", value=T, msg$get.types())
		if (length(error.type) > 0) .tag = ifelse(error.type[1] != "data.error", gsub(".", " ", error.type[1], fixed=T), "FAILURE")
	}
	if (is.null(.tag)) .tag = "WARNING"
	msg = msg$set(.prefix=paste0(toupper(.tag), ": "))$get.message()
	if (!globals$warning.message.state) msg = paste0("\n", msg)
	message(msg); globals$warning.message.state = T
}

## three general subclasses of errors:
## - fatal error: invalid program state, should not be reached regardless of user input; abort program, requires developer debugging
## - input error: invalid (user) input (eg. wrong file path or incorrect input file formatting); abort program, requires different input
## - data error: invalid data state, cannot proceed with current computation; program and data state remain valid, should generally be caught and handled by calling code
##   - printed as FAILURE when converted to warning
throw.error = function(...) {stop(create.message(..., .as.error=T))}

fatal.error = function(..., .sub.type=NULL) {throw.error(..., .sub.type=c(.sub.type, "fatal.error"))}
input.error = function(..., .sub.type=NULL) {throw.error(..., .sub.type=c(.sub.type, "input.error"))}
data.error = function(..., .sub.type=NULL) {throw.error(..., .sub.type=c(.sub.type, "data.error"))}

undefined.error = function(func.name, class.name) {fatal.error("function ", func.name, "() is not defined for class ", class.name, .sub.type="undefined.function")}

is.error = function(object, type=NULL) {return(has.type(object, c("error", type), R6=F))}
is.failure = function(object, type=NULL) {return(has.type(object, c("data.error", "error", type), R6=F))}

error.message = function(error, main.only=T) {return(ifelse(is.error(error), ifelse(has.type(error, "ConditionMessage"), ifelse(main.only, error$get.main(), error$get.message()), error$message), ""))}


## convert incoming errors of types listed in as.warning/as.error, as well as all data.error/error if all.failure/all.error are TRUE
## if error caught, prints error message and returns NULL; for error types in as.warning, print WARNING tag instead of original error tag
catchErrors = function(expr, as.warning=NULL, as.error=NULL, all.failures=F, all.errors=F) {
	catch.list = c(as.warning, as.error, if (all.failures) "data.error", if (all.errors) "error")
	catch.func = function(error) {
		caught = catch.list %in% class(error)
		if (any(caught)) throw.warning(error, .tag=if (any(catch.list[caught] %in% as.warning)) "WARNING")
		else throw.error(error)
		return(NULL)
	}
	tryCatch(expr, error=catch.func)
}

## execute expression, and return any errors as object instead of expression output if caught
## test result with is.error() for subsequent handling
# TODO: also capture LAVA generated warnings
muffleErrors = function(expr, allow.warnings=F) {
	if (!allow.warnings) return(suppressWarnings(muffleErrors(expr, allow.warnings=T)))
	else return(tryCatch(expr, error=function(err) {return(err)}))
}


create.message = function(..., .sub.type=NULL, .as.error=F) {args = list(...)
	is.msg = unlist(lapply(args, has.type, "ConditionMessage"))
	is.error = unlist(lapply(args, is.error));

	if (!any(is.msg)) {
		message = ConditionMessage$new(...)
		if (any(is.error)) {.sub.type = c("[r].error", .sub.type); .as.error = T}
	} else message = args[[which(is.msg)[1]]]$copy()$set(...)

	if (.as.error) .sub.type = c("simpleError", "error", "condition", .sub.type)
	message$set(.types=.sub.type)

	class(message) = unique(c(class(message)[1:which(class(message) == "R6")], message$get.types()))
	return(message)
}

## message consists of [PREFIX] [MAIN] ([LABELS]), followed by indented [LINES]
## all except [MAIN] are optional, [LABELS] are printed as semi-colon separated NAME: VALUE pairs (name is optional)
## input to [MAIN] is concatenated together, with named vectors in [MAIN] with name starting with "items" are pasted into a comma-separated, single-quoted list
## [MAIN] also pluralized tokens surrounded by %'s, see pluralize() function below
## additional data and types can be embedded in the message as well
ConditionMessage = R6::R6Class("ConditionMessage",
	cloneable = FALSE,
	private = list(
		main = "", prefix = "", labels = list(), lines = c(),
		types = c(), data = list(),

		parse.main = function(args) {
			if (length(args) > 0) {
				args = lapply(args, function(a) {if ("message" %in% names(a)) ifelse(has.type(a, "ConditionMessage"), a$get.main(), a$message) else a})
				private$main = pluralize(args)
			} else private$main = ""
		},

		parse.labels = function(labels) {
			if (!is.null(labels)) {
				if (length(labels) == 1 && is.na(labels)) labels = list()
				if (!is.list(labels)) labels = as.list(labels)
				private$labels = labels[!sapply(labels, is.null)]
			}
		},

		set.message = function() {
			if (length(private$labels) > 0) {
				parts = unlist(private$labels)
				if (!is.null(names(private$labels))) {
					named = names(private$labels) != ""
					parts[named] = paste0(names(private$labels)[named], ": ", parts[named])
				}
				label.str = paste0(" (", paste(parts, collapse="; "), ")")
			} else label.str = ""
			lines.str = if (length(private$lines) > 0) paste0("\n\t", private$lines, collapse="")

			self$message = paste0(private$prefix, private$main, label.str, lines.str, "\n")
		}
	),
	public = list(
		message = "", call = NULL,
		initialize = function(...) {self$set(...)},

		set = function(..., .prefix=NULL, .label=NULL, .lines=NULL, .types=NULL, .data=NULL) {
			if (length(list(...)) > 0) private$parse.main(list(...))
			if (!is.null(.label)) private$parse.labels(.label)
			if (!is.null(.prefix)) private$prefix = .prefix
			if (!is.null(.lines)) private$lines = unlist(.lines)

			if (!is.null(.types)) {
				private$types = unique(c(private$types, .types))
				error.index = grep(".*\\.error", private$types)
				if (length(error.index) > 1) private$types = private$types[-rev(error.index)[-1]]  ## remove competing error types (ie. any except last of form "error.[TYPE]")
			}
			if (!is.null(.data)) private$data = c(private$data, .data)

			private$set.message()
			invisible(self)
		},

		print = function(...) {cat(self$message)},

		get.main = function() {return(private$main)},
		get.message = function(newline=F) {return(ifelse(!newline, gsub("\n$", "", self$message), self$message))},

		get.data = function() {return(private$data)},
		get.types = function() {return(private$types)},
		has.type = function(type) {return(type %in% self$get.types())},

		copy = function() {
			out = ConditionMessage$new()
			for (main in c("self", "private")) {
				curr = self$.__enclos_env__[[main]]
				for (arg in names(curr)) {
					if (!(is.function(curr[[arg]]) || is.environment(curr[[arg]]))) out$.__enclos_env__[[main]][[arg]] = curr[[arg]]
				}
			}
			class(out) = class(self)
			return(out)
		}
	)
)


## concatenates input into single string, pluralizing tokens surrounded by %'s
## in general, messages should be written in singular form, eg. list("the %variable% ", items.and=variable.names, " %is% invalid or %does% not exist")
## - this will print as "the variable 'V1' is invalid or does not exist" if variable.names="V1"
## - or as "the variables 'V1' and 'V2' are invalid or do not exist" if variable.names=c("V1", "V2")
## NB: leading and trailing spaces (but not tabs) are removed, and multiple consecutive spaces are reduced to one

## special arguments:
## - .plural: determines if string should be singular or plural
##   - if not NULL, can be a TRUE/FALSE value, a number (set to TRUE if > 1), or a vector (set to TRUE if length > 1)
##   - if NULL, will be set to TRUE if any of the following special tokens are present
##     - items.count element of value > 1
##     - items.[WORD] element with length > 1
##     - any object with .plural attribute set to TRUE
## - items / items.[WORD]: vector, concatenated into comma-separated list of single-quoted items; last two items are separated by [WORD], if provided
## - items.count: printed directly (use .plural or object with .plural attribute to pluralize message without additional printing)

## pluralizing tokens:
## - %token%: no whitespace; pluralized by appending an s, or by lookup if token is listed in globals$plural.form in globals.R
## - %[token]%: same as %token%, except allows whitespace
## - %[+token]%: print token if plural, nothing if singular
## - %[-token]%: print token if singular, nothing if plural
## - %token[+part]% / %[+part]token%: will print tokenpart/parttoken if plural, token if singular
## - %token[-part]% / %[-part]token%: will print tokenpart/parttoken if singular, token if plural
## - %[?token_S/token_P]%: print token_S if singular, token_P if plural
pluralize = function(...) {args = flatten.arglist(..., filter.unnamed=F)
	plurify = function(word, as.plural) {ifelse(as.plural, ifelse(word %in% names(globals$plural.form), globals$plural.form[[word]], paste0(word, "s")), word)}

	## process all items.* tokens into single string, set make.plural value
	items = substr(names(args), 1, 5) == "items"
	if (any(items)) {
		make.plural = any(sapply(args[items], length) > 1) || (!is.null(args$items.count) && any(unlist(args[names(args) == "items.count"]) > 1))
		for (i in which(items & names(args) != "items.count"))	args[[i]] = paste0(" ", quote.vector(args[[i]], last=gsub("items\\.?", "", names(args)[i])), " ")
	} else make.plural = F

	## check for .plural argument (overrides previous value if present), and for any tokens with .plural attribute
	if (!is.null(names(args)) && !is.null(args$.plural)) {  ## override items counts
		make.plural = ifelse(!is.logical(args$.plural), length(args$.plural) > 1 || (is.numeric(args$.plural) && args$.plural > 1), args$.plural)
		args = args[names(args) != ".plural"]
	} else if (!make.plural) make.plural = any(sapply(args, function(x) {!is.null(attr(x, ".plural")) && attr(x, ".plural")}))

	## collapse into single string
	msg = paste0(unlist(args), collapse="")

	## parse %[ ]% special tokens
	elem = unique(stringr::str_extract_all(msg, "%[^%]*\\[[^%]+\\][^%]*%")[[1]])  ## %pref[@token]suf% component, pref/suf can be empty, with @ denoting the special character +/-/? denoting the mode
	for (e in elem) {
		token = gsub("^%\\S*(\\[[^%]+\\])\\S*%$", "\\1", e)   ## extract [@token] part
		token.mode = gsub("^\\[(.).*\\]$", "\\1", token)
		token.body = gsub("^\\[.(.*)\\]$", "\\1", token)

		if (token.mode == "?") {
			parts = unlist(strsplit(token.body, "/"))
			value = ifelse(make.plural, tail(parts, 1), head(parts, 1)) ## if no / in token, will equal token.body; if multiple slashes, uses string before first / or string after last /, disregards middle
		} else if (token.mode == "+") value = ifelse(make.plural, token.body, "")
		else if (token.mode == "-") value = ifelse(make.plural, "", token.body)
		else value = plurify(paste0(token.mode, token.body), make.plural)  ## %[token]% mode

		value = gsub("%", "", gsub(token, value, e, fixed=T))  ## replace [@token] in element and strip off %'s
		msg = gsub(e, value, msg, fixed=T)  ## replace %pref[@token]suf% component in message string
	}

	## pluralize regular %token% components; do lookup in global$plural.form list for 'token' to pluralize, append 's' if not present
	elem = unique(stringr::str_extract_all(msg, "%\\S+%")[[1]])
	for (e in elem) msg = gsub(e, plurify(substr(e, 2, nchar(e)-1), make.plural), msg, fixed=T)

	## strip leading/trailing spaced, and collapse consecutive spaces into one
	return(gsub("^ +", "", gsub(" +$", "", gsub(" +", " ", msg))))
}


## adds a .plural attribute; note that this is stripped by c() function
add.plural = function(object, is.plural) {return(add.info(object, .plural=is.plural))}




##############

## collapse input into string, quoting individual elements
quote.vector = function(..., quote="'", sep=", ", last="") {
	args = paste0(quote, c(...), quote)
	if (last != "" && length(args) > 1) args = c(head(args, -2), paste(tail(args, 2), collapse=paste0(" ", last, " ")))
	return(paste(args, collapse=sep))
}

quote.items = function(label, ..., quote="'", sep=", ") {
	if (length(c(...)) > 1) label = paste0(label, "s")
	return(paste0(label, " ", quote.vector(..., quote=quote, sep=sep)))
}

parentheses = function(value) {return(if (!is.null(value)) paste0("(", value, ")"))}


#####################




## wrapper class for easier printing of LAVA objects, build up in parts then convert to string to print
## add.line() pastes arbitrary content together into single line
## add.parameter() prints parameter as "name: value", prints nothing if value is NULL
## add.list() prints "name: length" followed by "parameter: value" pairs on separate lines add indent+1 from the named list used as input
## - nothing is printed if the list is empty unless show.empty=T
## - if values.only=T or 'values' is not a named list, only the values are printed on each line
## parameter.list() expects arbitrary number of 'parameter=values' arguments, printed as "parameter: values[1], values[2], ..."
ObjectPrinter = R6::R6Class("ObjectPrinter",
	private = list(name = NULL, lines = NULL),
	public = list(
		initialize = function(...) {self$set.name(...)},

		set.name = function(...) {private$name = paste0(...); invisible(self)},
		add.line = function(..., add.space=T, indent=1) {
			prefix = paste0(rep("  ", times=indent), collapse="")
			private$lines = c(private$lines, paste0(prefix, paste(..., sep=ifelse(add.space, " ", ""))))
			invisible(self)
		},
		add.parameter = function(name, value, indent=1) {if (!is.null(value)) self$add.line(name, ": ", value, add.space=F, indent=indent); invisible(self)},
		add.list = function(name, values, values.only=F, show.size=T, indent=1) {
			if (length(values) > 0 || show.size) self$add.line(name, if (show.size) paste0(": ", length(values)), add.space=F, indent=indent)
			if (values.only || !is.list(values) || is.null(names(values))) {
				values = as.list(values)
				for (i in seq_along(values)) self$add.line(values[[i]], indent=indent+1)
			} else self$parameter.list(values, indent=indent+1)
			invisible(self)
		},
		parameter.list = function(..., indent=1) {
			args = flatten.arglist(...); args = args[!sapply(args, is.null)]
			for (arg in names(args)) self$add.line(arg, ": ", paste0(args[[arg]], collapse=", "), add.space=F, indent=indent)
			invisible(self)
		},

		to.string = function() {paste0("LAVA ", private$name, " object\n", paste0(private$lines, collapse="\n"), "\n")}
	)
)

