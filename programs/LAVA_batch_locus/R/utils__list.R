
## if multiple values provided: creates a list of all the input items
## if no additional values provided, and 'values' is not already a list:
## - if length(names) == 1, creates a list with 'values' as it's first element
## - otherwise, creates a list with each element of 'values' as a separate list item
## if .drop.null=T, NULL elements are removed after constructing the list, ie. length of names should match total number of input values provided
named.list = function(names, values, ..., .drop.null=T) {args = list(...)
	if (is.null(names)) names = character(0)
	if (length(args) > 0) values = list(values, ...)
	else if (!is.list(values)) values = if (length(names) == 1) list(values) else as.list(values)
	if (length(names) != length(values)) input.error("number of names does not match number of values in function named.list()")
	names(values) = names
	return(if (.drop.null && length(values) > 0) values[!unlist(lapply(values, is.null))] else values)
}


## clean named list: remove duplicates (retain first/last), and remove unnamed and NULL elements
clean.list = function(input, duplicates=c("last", "first"), drop.null=T) {
	if (!is.null(input) && !is.list(input)) input.error("input to clean.list() is not a list")

	if (length(input) > 0 && !is.null(names(input))) {
		input = input[names(input) !=	""]
		input = input[!duplicated(names(input), fromLast=match.arg(duplicates) == "last")]
		if (drop.null && length(input) > 0) input = input[!sapply(input, is.null)]
	}
	if (length(input) == 0 || is.null(names(input))) {input = list(); attr(input, "names") = character(0)}
	return(input)
}


## remove zero-length elements from list
trim.list = function(lst) {return(lst[unlist(lapply(lst, length)) > 0])}

## append input in ... as additional list elements (rather than appending content like c() / append())
list.append = function(lst, ...) {add = list(...)
	if (length(add) > 0) lst[1:length(add) + length(lst)] = add
	return(lst)
}

## from a list of lists, and element names in ... to extact, will create a named (element names) list with for each element a list containing the corresponding contents of the inner list
##   eg. if lst = list(x1 = list(names=[names1], values=[values1]), x2 = list(names=[names2], values=[values2])), then will output:
##   list(names = list(x1=[names1], x2=[names2]), values = list(x1=[values1], x2=[values2]))
## element names can contain $ signs, eg. "values$rho" will extract lst[[i]]$values$rho for every index i; output element name will be the last part only, ie. 'rho' in the example
## names of the input list 'lst' (ie. x1, x2 above) are preserved if present
## if trim=F, each inner list will be the same length as 'lst', if trim=T empty elements are removed
## if flatten=T, each element will be flattened to a vector using unlist()
## if collapse=T and only a single element was requested, the corresponding internal list will be returned
##   eg. list(x1=[values1], x2=[values2]) instead of list(values = list(x1=[values1], x2=[values2]))
list.extract = function(lst, ..., trim=F, flatten=F, collapse=T) {
	output = list();
	for (label in c(...)) output[[label]] = lapply(lst, function(l) {eval(parse(text=paste0("l$", label)))})
	names(output) = gsub("^.*\\$([^\\$]*)$", "\\1", names(output))
	if (trim == "auto") trim = !is.null(names(lst))
	if (trim) output = lapply(output, trim.list)
	if (flatten) output = lapply(output, unlist)
	if (collapse && length(output) == 1) return(output[[1]])
	else return(output)
}

## merge named lists and additional named arguments, creating a single named list
## named arguments are treated as a list, ie. list.merge(a, b=1) is equivalent to list.merge(a, list(b=1))
## unnamed non-list arguments are discarded, as are named NULL elements unless drop.null=F
## for entries with the same name, the last value is used
list.merge = function(..., drop.null=T) {
	args = list(...)
	if (length(args) > 0) {
		if (!is.null(names(args))) {for (i in which(names(args) != "")) args[[i]] = add.names(list(args[[i]]), names=names(args)[i])}
		args = lapply(args[sapply(args, is.list)], clean.list, drop.null=drop.null)
	}

	if (length(args) > 1) {
		out = args[[1]]
		for (i in 2:length(args)) out[names(args[[i]])] = args[[i]]
		return(out)
	} else return(if (length(args) == 1) args[[1]] else clean.list(list()))
}


