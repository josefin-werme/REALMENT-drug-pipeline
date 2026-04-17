## sum columns of input matrix (considerably faster than apply(data, 1, sum) for large inputs)
column.sum = function(data) {return((data %*% rep(1, ncol(data)))[,1])}

## apply logical function to boolean matrix (considerably faster than apply(bool.data, 1, op) for large inputs)
column.none = function(bool.data) {return(column.sum(bool.data) == 0)}
column.any = function(bool.data) {return(column.sum(bool.data) > 0)}
column.all = function(bool.data) {return(column.sum(bool.data) == ncol(bool.data))}






## subset matrix by row/colnames; cols/rows values of NULL return all rows/columns
## maintains row/colnames in output and correct dimensions even if size of one or both dimensions reduced to zero
matrix.subset = function(M, cols=NULL, rows=NULL) {
	if (!is.null(rows)) {
		if (is.null(rownames(M))) input.error("in matrix.subset(), cannot select rows if rownames is NULL")
		if (!all(rows %in% rownames(M))) input.error("in matrix.subset(), unknown rows selected")
		row.index = match(rows, rownames(M))
	} else row.index = 1:nrow(M)
	if (!is.null(cols)) {
		if (is.null(colnames(M))) input.error("in matrix.subset(), cannot select cols if colnames is NULL")
		if (!all(cols %in% colnames(M))) input.error("in matrix.subset(), unknown cols selected")
		col.index = match(cols, colnames(M))
	} else col.index = 1:ncol(M)

	return(matrix(M[row.index,col.index], nrow=length(row.index), ncol=length(col.index), dimnames=list(rownames(M)[row.index], colnames(M)[col.index])))
}






##############################
###   COVARIANCE MATRICES  ###
##############################


## create diagonal matrix, optionally set row/column names
## argument 'diagonal' is always treated as a vector, even if length is one
diag.matrix = function(diagonal, add.names=NULL) {
	M = diag(diagonal, nrow=length(diagonal))
	if (!is.null(add.names)) {
		if (length(add.names) != length(diagonal)) input.error("add.names argument for diag.matrix() does not match size of specified matrix")
		colnames(M) = rownames(M) = add.names
	}
	return(M)
}


## elements should contain correlations in lower diagonak, in column order; multiple input elements / vectors are appended together
correlation.matrix = function(...) {
	values = c(...)
	if (length(values) > 0) {
		size = (sqrt(8*length(values) + 1) - 1) / 2
		if (size == round(size)) {
			M = diag(round(size)+1)/2; M[lower.tri(M)] = values; M = M + t(M)
			if (any(eigen(M)$values < 0)) input.error("correlation matrix is not positive-semidefinite")
			return(M)
		} else input.error("incorrect number of correlations specified")
	} else return(matrix(1, 1, 1))
}


## create correlation matrix corresponding to model with Y = X*beta + err, with given r2
regression.matrix = function(cov.x, beta, r2) {
	cov.xy = cov.x %*% beta
	var.y = ifelse(r2 > 0, t(beta) %*% cov.xy / r2, 1)
	return(cov2cor(rbind(cbind(cov.x, cov.xy), c(cov.xy, var.y))))
}



## mirror square matrix values and names across the diagonal, and optionally set new row/column names
make.symmetric = function(M, add.names=NULL, use=c("lower", "upper")) {
	if (is.null(dim(M)) || nrow(M) != ncol(M)) input.error("matrix input for make.symmetric() is not square")
	if (match.arg(use) == "upper") M = t(M)

	if (!is.null(add.names)) {
		if (length(add.names) != nrow(M)) input.error("add.names argument for make.symmetric() does not match size of matrix input")
		rownames(M) = add.names
	}

	colnames(M) = rownames(M)
	M[upper.tri(M)] = t(M)[upper.tri(M)]
	return(M)
}


## wrapper to retrieve underlying matrix if CovarianceMatrix object
matrix.contents = function(M) {return(if (has.type(M, "CovarianceMatrix")) M$get() else M)}


## convert input matrix to CovarianceMatrix (if not already); sets names to numeric string dummies if NULL
as.covariance = function(M, names=NULL, type=NULL) {
	if (!has.type(M, "CovarianceMatrix") && is.null(colnames(M)) && is.null(rownames(M)) && is.null(names) && !is.null(dim(M)) && nrow(M) > 0) names = 1:nrow(M)
	return(check.covariance(M, names=names, type=type))
}


## create CovarianceMatrix object with given names, variance vector, and constant covariance value
create.covariance = function(names, covar.value=0, variances=1, type="covariance matrix") {
	size = length(names)
	if (any(duplicated(names))) input.error("cannot create ", type, ", argument 'names' contains duplicate entries")
	if (length(variances) > 1 && length(variances) != size) input.error("cannot create ", type, ", number of names does not match number of variance values")

	M = matrix(covar.value, size, size)
	if (size > 0) diag(M) = variances
	CovarianceMatrix$new(M, names=names, type=type)
}


## check if covariance matrix is valid; matrix must be named unless explicit 'names' argument is provided; returns CovarianceMatrix object
## allow.NA=T setting only allows NA covariance values, not variances
check.covariance = function(M, names=NULL, allow.NA=T, type=NULL) {
	covar = CovarianceMatrix$new(M, names=names, type=type)
	if (!allow.NA && any(is.na(covar$get()))) input.error(ifelse(is.null(type), "covariance matrix", type), " contains missing values")
	return(covar)
}



## if 'names' argument is set, it will overwrite any existing col/row names

## STATUS: stable, shallow
CovarianceMatrix = R6::R6Class("CovarianceMatrix",
	private = list(
		data = NULL, type = "covariance matrix", tolerance = 1e-10,

		validate = function(M, names) {
			if (is.null(M) || !is.matrix(M)) input.error(private$type, " is not a matrix")
			if (is.null(dim(M)) || nrow(M) != ncol(M)) input.error(private$type, " is not a square matrix")

			if (nrow(M) > 0) {
				if (!isSymmetric(M, tol=private$tolerance, check.attributes=F)) input.error(private$type, " is not symmetrical")
				if (any(is.na(diag(M)) | diag(M) <= 0)) input.error(private$type, " contains missing, negative or zero variance values")
				if (any(abs(cov2cor(M)) > 1, na.rm=T)) input.error(private$type, " contains correlations exceeding 1 or -1")

				if (is.null(names)) {
					if (!is.null(colnames)) {
						if (!is.null(rownames(M)) && !all(colnames(M) == rownames(M))) input.error("row and column names of ", private$type, " do not match")
						names = colnames(M)
					} else names = rownames(M)
					if (is.null(names)) input.error("no row/column names provided for ", private$type)
				}
				if (length(names) != nrow(M)) input.error("length of 'names' argument does not match dimension of ", private$type)
				if (any(duplicated(names))) input.error(private$type, " contains duplicate row/column names")
				return(add.dimnames(M, names))
			} else return(matrix(0,0,0))
		},

		modified.copy = function(M, names=NULL) {
			if (is.null(names) && is.null(colnames(M)) && is.null(rownames(M)) && nrow(M) > 0 && ncol(M) > 0) {
				if (nrow(M) != self$size()) fatal.error("no valid names available in in modified.copy()")
				names = self$get.names()
			}
			return(CovarianceMatrix$new(M, names=names, type=private$type, tolerance=private$tolerance))
		}
	),
	public = list(
		initialize = function(M, names=NULL, type=NULL, tolerance=NULL) {
			if (has.type(M, "CovarianceMatrix")) {
				if (is.null(type)) type = M$get.type()
				if (is.null(tolerance)) tolerance = M$.__enclos_env__$private$tolerance
				M = M$get()
			}
			if (!is.null(type)) private$type = type
			if (!is.null(tolerance)) private$tolerance = tolerance
			private$data = private$validate(M, names=names)
		},

		print = function(...) {
			cat("LAVA CovarianceMatrix object", if (private$type != "covariance matrix") parentheses(private$type), "\n\n")
			print(private$data)
		},

		size = function() {return(nrow(private$data))},
		get = function() {return(private$data)},

		get.names = function() {return(if (nrow(private$data) > 0) rownames(private$data) else character(0))},
		get.type = function() {return(private$type)},

		standardize = function() {
			if (any(abs(diag(private$data) - 1) > private$tolerance)) {
				M = cov2cor(private$data); diag(M) = 1
				return(private$modified.copy(M))
			} else (return(self))
		},

		scale = function(variances, standardize=F) {
			if (length(variances) != self$size()) input.error("cannot rescale covariance matrix, incorrect number of variances provided")
			if (any(is.na(variances) | variances <= 0)) input.error("cannot rescale covariance matrix, variance values are invalid")
			M = if (standardize) cov2cor(private$data) else private$data
			S = diag(sqrt(variances), ncol=length(variances))
			return(private$modified.copy(S %*% M %*% S))
		},

		## also places elements in input order
		subset = function(names, ignore.missing=F) {
			missing = !(names %in% self$get.names())
			if (any(missing)) {
				if (ignore.missing) names = names[!missing]
				else input.error(private$type, " does not contain %element% ", items.and=names[missing])
			}
			if (self$size() > 0) return(private$modified.copy(private$data[names,names,drop=F]))
			else return(self)
		},

		## expand with additional rows/columns, with constant covariance value
		expand = function(names, covar.value=0, variances=1, ignore.existing=F) {
			if (any(duplicated(names))) input.error("cannot expand ", private$type, ", argument 'names' contains duplicate entries")
			if (length(names) == 0) input.error("cannot expand ", private$type, ", argument 'names' is empty")
			if (length(variances) > 1 && length(names) != length(variances)) input.error("cannot expand ", private$type, ", number of names does not match number of variance values")

			existing = names %in% self$get.names()
			if (any(existing)) {
				if (ignore.existing) {
					if (all(existing)) return(self)
					names = names[!existing]
					if (length(variances) > 1) variances = variances[!existing]
				} else input.error("cannot expand ", private$type, ", %element% ", items.and=names[existing], " already %exist[-s]%")
			}

			total = self$size() + length(names); add.range = self$size() + 1:length(names)
			M = add.dimnames(matrix(covar.value, total, total), c(self$get.names(), names))
			diag(M)[add.range] = variances
			M[-add.range,-add.range] = private$data
			return(private$modified.copy(M))
		},

		## return copy with new type
		set.type = function(type) {out = self$clone(); out$.__enclos_env__$type = type; return(out)},

		## return copy with new names
		set.names = function(names) {return(private$modified.copy(private$data, names))},

		## return copy with new values; if no names are set on M, will use existing names
		set.values = function(M, names=NULL) {
			if (is.null(names) && is.null(colnames(M)) && is.null(rownames(M))) names = colnames(private$data)
			out = self$clone()
			out$.__enclos_env__$private$data = private$validate(M, names=names)
			return(out)
		}
	)
)

