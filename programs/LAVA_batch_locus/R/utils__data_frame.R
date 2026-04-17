
## merge list of data frames, aligning column order and filling in missing columns
fill.rowmerge = function(...) {data = list(...)
	if (length(data) == 1 && is.list(data[[1]]) && !is.data.frame(data[[1]])) data = data[[1]]
	data = data[!sapply(data, is.null)]
	columns = consensus.order(lapply(data, names))
	if (!all(sapply(data, ncol) == length(columns)) || !all(sapply(data, function(d) {all(names(d) == columns)}))) {
		for (i in seq_along(data)) {
			missing = !(columns %in% names(data[[i]]))
			if (any(missing)) data[[i]][columns[missing]] = NA
			data[[i]] = data[[i]][columns]
		}
	}
	return(do.call(rbind, data))
}


## select rows in data with specific combination of values
## 'values' can either be a vector of length equal to the number of columns in 'data', or a named list, with names corresponding to column names of 'data'
## if multiple values are specified for a column, column will match if it equals any of those values
column.equals = function(data, values) {
	if (!is.list(values)) {
		if (length(values) != ncol(data)) input.error("number of values must match number of columns, in column.equals()")
		values = add.names(as.list(values), names(data))
	}
	if (is.null(names(values)) || !all(names(values) %in% names(data))) input.error("invalid value names, in column.equals()")

	keep = rep(T, nrow(data))
	for (var in names(values)) {
		if (length(values[[var]]) > 1) keep = keep & data[[var]] %in% values[[var]]
		else keep = keep & data[[var]] == values[[var]]
	}
	return(keep)
}



## paste columns of input data.frame together
paste.columns = function(data, columns, sep=" ") {
	if (length(columns) > 1) {
		if (is.numeric(columns)) columns = names(data)[columns]
		expr.str = paste0("paste(", paste0("data[[\"", columns, "\"]]", collapse=", "), ", sep=\"", sep, "\")")
		return(eval(parse(text=expr.str)))
	} else return(data[[columns]])
}


