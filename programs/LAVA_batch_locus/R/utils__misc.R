## add names to object
add.names = function(obj, names, check.matrix=T) {
	if (check.matrix && is.matrix(obj)) return(add.dimnames(obj, row.names=names))
	names(obj) = names; return(obj)
}
add.colnames = function(obj, names) {colnames(obj) = names; return(obj)}
add.rownames = function(obj, names) {rownames(obj) = names; return(obj)}
add.dimnames = function(obj, row.names, col.names=NULL) {
	rownames(obj) = row.names
	colnames(obj) = if (is.null(col.names)) row.names else col.names
	return(obj)
}


## update names specified in ... as old.name=new.name value pairs, if name is present on object (ignored otherwise)
## ignores case, does not verify if name already exists
change.names = function(obj, ...) {
	update = flatten.arglist(...); names(update) = tolower(names(update))
	for (old in names(update)) names(obj)[tolower(names(obj)) == old] = update[[old]]
	return(obj)
}



## set named values in ... as attributes on the input object
add.info = function(object, ...) {
	args = flatten.arglist(...)
	for (var in names(args)) attr(object, var) = args[[var]]
	return(object)
}


## filter out NA entries in a list or vector, returns untouched if neither; also filters NULL values from lists if drop.null=T
filter.NA = function(values, drop.null=T) {
	if (length(values) > 0) {
		if (is.list(values)) return(values[!sapply(values, function(v) {(drop.null && is.null(v)) || (length(v) == 1 && is.na(v))})])
		else if (is.vector(values)) return(values[!is.na(values)])
	}
	return(values)
}


## tries to create a consensus ordering of vectors of input values, returns an ordered vector of all unique values in the input
## primarily intended for use with vectors that have the same ordering but have different values missing
consensus.order = function(values) {
	range = unique(unlist(values))
	if (!all(sapply(values, length) == length(range)) || !all(sapply(values, function(v) {all(v == range)}))) {
		range = c()
		while (length(values) > 1) {
			current = unlist(lapply(values, head, 1))
			if (length(unique(current)) > 1) {
				first.only = sapply(current, function(c) {all(sapply(values, function(vals) {!(c %in% vals[-1])}))})
				if (any(first.only)) current = current[first.only]
				if (length(unique(current)) > 1) {
					tbl = table(current)
					current = current[current %in% names(tbl)[tbl == max(tbl)]]
				}
			}
			range = c(range, current[1])
			values = trim.list(lapply(values, function(vals) {vals[vals != current[1]]}))
		}
		if (length(values) == 1) range = c(range, values[[1]])
	}
	return(range)
}


## create list of unique pairs of values from two input vectors
## if ignore.order=T, pairs are considered duplicates if they contain the same values regardless of order
## if allow.equal=F, pairs are removed if their values are the same
make.pairs = function(input1, input2=NULL, allow.equal=F, ignore.order=T, sort=F) {
	if (is.null(input2)) input2 = input1
	if (is.list(input1) || is.list(input2)) fatal.error("input to function make.pairs() must be vectors, not lists")
	if (length(input1) > 0 && length(input2) > 0) {
		pairs = unlist(lapply(input1, function(v1) {lapply(input2, function(v2) {c(v1,v2)})}), recursive=F)
		if (sort) pairs = lapply(pairs, sort)
		if (!allow.equal) pairs = pairs[!sapply(pairs, function(p) {p[1] == p[2]})]
		if (length(pairs) > 0) pairs = pairs[!duplicated(if (ignore.order) lapply(pairs, sort) else pairs)]
		return(pairs)
	} else return(list())
}


## create a list of all possible subsets (ignoring order) of the vector of input values
## if allow.duplicates=T, will consider duplicate values as distinct; otherwise, will remove them before creating subsets
make.subsets = function(values, allow.duplicates=F) {
	if (!allow.duplicates) values = unique(values)
	if (length(values) > 1) {
		size = length(values)
		index = as.list(1:size); current = index
		for (s in 2:size) {
			current = unlist(lapply(current, function(ind) {m = max(ind); if (m < size) lapply((m+1):size, function(x) {c(ind, x)})}), recursive = F)
			index = c(index, current)
		}
		values = sort(values)
		return(lapply(index, function(ind) {values[ind]}))
	} else return(if (length(values) == 1) list(values) else list())
}


## create data.frame with all combinations of named input vectors or input named list
combinations = function(...) {
	args = flatten.arglist(...); values = NULL
	for (param in names(args)) {
		if (!is.null(values)) values = cbind(values[rep(1:nrow(values), each=length(args[[param]])),,drop=F], args[[param]])
		else values = data.frame(args[[param]])
	}
	names(values) = names(args); rownames(values) = NULL
	return(values)
}


## for vector of strings, match against pattern containing tokens with format {TOKEN_NAME}, ignoring regular expression characters in pattern string
## return data.frame with matching strings in first column, and values for each token in subsequent named columns
## if same token occurs multiple times in input pattern, only strings where token values are all the same for that token are retained
## token names are returned in lower case; if 'allowed' is set, only tokens with specified names are extracted (if present; case insensitive)
extract.tokens = function(pattern, values, allowed=NULL) {
	tokens = stringr::str_extract_all(pattern, "\\{[^\\}]+\\}")[[1]]
	if (length(tokens) > 0) {
		token.names = tolower(gsub("[\\{|\\}]", "", tokens))
		chunks = stringr::str_escape(stringr::str_split(pattern, "\\{[^\\}]+\\}")[[1]])
		if (length(chunks) != length(tokens) + 1) fatal.error("unable to process pattern '", pattern, "' in function extract.tokens()")

		keep = if (!is.null(allowed)) token.names %in% allowed else rep(T, length(tokens))
		insert = ifelse(keep, "(.+)", stringr::str_escape(tokens))
		token.names = token.names[keep]
		match.str = paste0(chunks[1], paste(paste0(insert, chunks[-1]), collapse=""))

		matched = data.frame(stringr::str_match(values, match.str)); names(matched) = c("", token.names)
		matched = matched[!apply(is.na(matched), 1, any),,drop=F]

		for (dup in unique(token.names[duplicated(token.names)])) matched = matched[apply(matched[,names(matched) == dup], 1, function(v) {length(unique(v))}) == 1,]
		matched = matched[,!duplicated(names(matched)),drop=F]
	} else matched = data.frame(input = values[values == pattern])
	return(matched)
}








