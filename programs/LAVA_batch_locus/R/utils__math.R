#' @include utils__input_validation.R
#' @include utils__misc.R
#' @include utils__matrix.R


## check invertibility of symmetric matrix M
is.invertible = function(M, threshold=globals$decomposition.eigen.threshold) {!is.error(muffleErrors(eigen.invert(M, threshold=threshold)))}


## invert symmetric matrix M using eigendecomposition, checking invertibility (NB: symmetry of M is assumed)
eigen.invert = function(M, threshold=globals$decomposition.eigen.threshold) {
	M = matrix.contents(M)
	decomp = eigen(M); no.vars = ifelse(!is.null(dim(M)), nrow(M), 1)
	if (any(decomp$values / sum(decomp$values) < threshold / no.vars)) data.error("matrix is not invertible")
	if (no.vars > 1) return(decomp$vectors %*% diag(1/decomp$values) %*% t(decomp$vectors))
	else return(matrix(1/M, 1, 1))
}

## compute root R of positive-semidefinite matrix M, such that M = R %*% t(R)
## does not check symmetry or positive-semidefinite status, truncates any negative eigenvalues to zero
matrix.root = function(M) {
	decomp = eigen(matrix.contents(M))
	decomp$values[decomp$values < 0] = 0
	return(decomp$vectors %*% diag(sqrt(decomp$values)))
}


## fit regression model from covariance matrix M, with outcome and predictor arguments specified in numeric indices
## by default, uses all other variables as predictors
fit.regression = function(M, outcome, predictors=NULL) {
	M = as.covariance(M)$get()
	if (is.null(predictors)) predictors = which(1:nrow(M) != outcome)
	if (outcome %in% predictors || length(predictors) == 0 || !all(c(outcome, predictors) %in% 1:nrow(M))) input.error("invalid input for function fit.regression")
	My = M[outcome,outcome]; Mxy = M[predictors,outcome,drop=F]
	Mx.inv = muffleErrors(eigen.invert(M[predictors,predictors,drop=F]))
	if (is.error(Mx.inv)) data.error("predictor covariance matrix is not invertible")

	out = list(coefficients = Mx.inv %*% Mxy, residual = My - t(Mxy) %*% Mx.inv %*% Mxy)
	out$r2 = 1 - out$residual / My
	return(out)
}



## build up linear structural model from input terms, defining every endogenous variable as a linear combination of other variables, and explained variance r2
## specified as formula terms, eg. y1 = x1 + 3*x2 - 5*x3
##   variable names are interpreted case insensitive, and have to start with a letter (can contain periods and numbers)
## variables are all defined as having a variance of 1, ie. weights in formula are -relative-
## any variables not defined in model are assumed to be exogenous, and -INDEPENDENT- of all other exogenous variables
StructuralModel = R6::R6Class("StructuralModel",
	private = list(
		settings = list( ## ObjectSettings object
			tolerance = add.info(1e-5, type="numeric", range=c(0,1))  ## relative variance threshold to set elements to 0
		),
		input.specification = list(),  ## named list of variable specifications, with value: list(name, r2, predictors=data.frame(variable, scale))
		model.data = list(
			exogenous = NULL,   ## vector of exogenous variables
			equations = list(),   ## named list of endogenous variables, with equation specified in terms of exogenous variables (with explicit residual) rescaled to a variance of one: data.frame(variable, scale)
			correlations = NULL,  ## CovarianceMatrix object, with correlation structure of all variables in model
			selection = NULL   ## vector of variables to show in output (show all if NULL)
		),

		## check phenotype names against available and resolve partial/wildcards; returns NULL if input is empty
		check.phenotypes = function(...) {
			private$resolve.model()
			return(if (length(c(...)) > 0) phenotypes = check.phenotypes(c(...), c(private$model.data$exogenous, names(private$model.data$equations)), label="structural model"))
		},

		## checks phenotypes, and intersects with selection (if any); returns NULL if empty input and no selection set
		check.selection = function(...) {
			phenotypes = private$check.phenotypes(...)
			if (length(phenotypes) > 0 && length(private$model.data$selection) > 0) return(intersect(phenotypes, private$model.data$selection))
			else return(c(phenotypes, private$model.data$selection))
		},

		mark.updated = function() {private$model.data = list()},
		check.variable = function(var) {if (!grepl("^[a-z][a-z0-9\\.]*$", var)) input.error("invalid variable name '", var, "'") else return(var)},

		parse.token = function(token) {
			scale.str = gsub("[a-z].*$", "", token)
			if (nchar(scale.str) > 0) {
				if (scale.str == "-") scale = -1
				else scale = suppressWarnings(as.numeric(gsub("\\*$", "", scale.str)))
				if (is.na(scale)) input.error("invalid token '", token, "'")
			} else scale = 1
			return(list(variable = private$check.variable(substring(token, nchar(scale.str)+1)), scale = scale))
		},

		parse.formula = function(formula) {
			parts = unlist(strsplit(tolower(gsub("\"|\\s", "", formula)), "|", fixed=T))  ## strips out quote characters and whitespace, converts to lower case
			if (length(parts) != 2) input.error("formula contains %[?no/multiple]% '|' %symbol%", .plural=length(parts) > 2)
			if (parts[2] == "") input.error("formula contains no predictors")
			outcome = private$check.variable(parts[1])

			tokens = unlist(strsplit(gsub("\\++", "\\+", gsub("-", "+-", parts[2], fixed=T)), "+", fixed=T))
			predictors = do.call(rbind, lapply(tokens, function(token) {as.data.frame(private$parse.token(token))}))
			IF.ANY(duplicate=duplicated(predictors$variable), THEN=input.error("formula contains duplicate %predictor% ", items.and=unique(predictors$variable[duplicated])))
			if (outcome %in% predictors$variable) input.error("formula contains outcome as predictor")

			return(list(outcome = outcome, predictors = predictors[order(predictors$variable),]))
		},

		resolve.model = function() {
			if (length(private$input.specification) == 0 || length(private$model.data) > 0) return()

			specification = private$input.specification[order(names(private$input.specification))]
			variables = sort(unique(c(names(specification), unlist(lapply(specification, function(v) {v$predictors$variable})))))
			exogenous = variables[!(variables %in% names(specification))]

			endogenous = c(); remainder = variables[!(variables %in% exogenous)]
			while (length(remainder) > 0) {  ## order endogenous variables by dependency, check for circularity
				resolved = sapply(specification[remainder], function(spec, reference) {return(all(spec$predictors$variable %in% reference))}, reference=c(exogenous, endogenous))
				if (any(resolved)) {
					endogenous = c(endogenous, remainder[resolved])
					remainder = remainder[!resolved]
				} else input.error("model contains circularities, could not resolve")
			}

			equations = list()
			for (v in endogenous) {
				predictors = specification[[v]]$predictors
				if (specification[[v]]$r2 > 0) {
					## replace all endogenous components with corresponding equation in terms of independent exogenous variables
					replace = which(predictors$variable %in% endogenous)
					for (i in replace) {
						add = equations[[predictors$variable[i]]]
						add$scale = add$scale * predictors$scale[i]
						predictors = rbind(predictors, add)
					}
					predictors = predictors[!(predictors$variable %in% endogenous),]
					predictors = predictors[order(grepl("^error::", predictors$variable), predictors$variable),]

					## collapse multiples
					while (any(duplicated(predictors$variable))) {
						curr.index = which(predictors$variable == predictors$variable[duplicated(predictors$variable)][1])
						predictors$scale[curr.index[1]] = sum(predictors$scale[curr.index])
						predictors = predictors[!(1:nrow(predictors) %in% curr.index[-1]),]
					}
				} else predictors$scale = 0

				## rescale to variance of one, remove redundant, and add explicit residual
				if (all(predictors$scale^2 <= private$settings$get("tolerance"))) predictors = predictors[-(1:nrow(predictors)),]
				else predictors = predictors[predictors$scale^2 / sum(predictors$scale^2) > private$settings$get("tolerance"),]

				if (nrow(predictors) > 0) {
					scale = specification[[v]]$r2 / sum(predictors$scale^2)
					predictors$scale = predictors$scale * sqrt(scale)
				}

				residual = ifelse(nrow(predictors) > 0, 1 - specification[[v]]$r2, 1)
				if (residual > 0) {
					error = paste0("error::", v)
					predictors[nrow(predictors)+1,] = list(variable = error, scale = sqrt(residual))
					exogenous = c(exogenous, error)
				}

				equations[[v]] = predictors
			}
			exogenous = exogenous[order(grepl("^error::", exogenous), exogenous)]

			## set implied correlation matrix
			pattern = matrix(0, nrow=length(exogenous), ncol=length(equations))
			for (i in seq_along(equations)) pattern[match(equations[[i]]$variable, exogenous),i] = equations[[i]]$scale
			M = rbind(cbind(diag(length(exogenous)), pattern), cbind(t(pattern), t(pattern) %*% pattern))

			private$model.data = list(exogenous = exogenous, equations = equations, correlations = as.covariance(M, names=c(exogenous, names(equations))))
		}
	),
	public = list(
		initialize = function(...) {private$settings = validate.settings(self, private$settings, ...)},
		add = function(formula, r2=1) {
			if (is.null(match.call()$formula)) input.error("no formula provided")
			if (r2 < 0 || r2 > 1) input.error("r2 value must be between 0 and 1")
			formula = private$parse.formula(deparse(match.call()$formula))
			private$input.specification[[formula$outcome]] = list(name = formula$outcome, r2 = r2, predictors = formula$predictors)
			private$mark.updated()
			invisible(self)
		},

		remove = function(...) {variables[names(variables) %in% tolower(c(...))] = NULL; private$mark.updated()},

		no.pheno = function(mode=c("all", "exogenous", "endogenous", "selection")) {return(length(self$phenotypes(match.arg(mode))))},
		phenotypes = function(mode=c("all", "exogenous", "endogenous", "selection")) {
			private$resolve.model(); mode = match.arg(mode)
			if (mode == "selection" && length(private$model.data$selection) > 0) return(private$model.data$selection)
			else return(c(if (mode != "endogenous") private$model.data$exogenous, if (mode != "exogenous") names(private$model.data$equations)))
		},

		get.specification = function() {return(private$input.specification)},
		get.model = function(...) {
			phenotypes = private$check.selection(...)

			if (length(phenotypes) > 0) {
				model = private$model.data
				model$equations[!(names(model$equations) %in% phenotypes)] = NULL
				model$exogenous = model$exogenous[model$exogenous %in% c(phenotypes, unlist(list.extract(model$equations, "variable")))]
				model$correlations = model$correlations$subset(c(model$exogenous, names(model$equations)))
				return(model)
			} else return(private$model.data)
		},

		get.correlations = function(...) {
			phenotypes = private$check.selection(...)
			return(if (length(phenotypes) > 0) private$model.data$correlations$subset(phenotypes) else private$model.data$correlations)
		},

		## set selection of phenotypes to show in output; is reset when updating model
		set.selection = function(...) {private$model.data$selection = if (length(c(...)) > 0) sort(private$check.phenotypes(...)); invisible(self)}
	)
)


## generate draws of matrices with means of zero and specified dependence across columns (covariance, or inproduct of matrix with itself)
## if exact=T, the product or (sample) covariance of the draw will be exactly equal to the expected value
## if weights are provided, rows are multiplied by the corresponding weight
MatrixSampler = R6::R6Class("MatrixSampler",
	private = list(
		no.rows = NULL, weights = NULL,
		settings = list(exact = F, mode = "covariance"),
		expectation.root = NULL
	),
	public = list(
		initialize = function(expectation, nrow, weights=NULL, exact=F, mode=c("covariance", "product")) {
			if (!is.null(weights)  && length(weights) != nrow) input.error("number of weights does not match number of rows")
			private$no.rows = nrow; private$weights = weights
			private$settings = list(exact = exact, mode = match.arg(mode))
			private$expectation.root = t(matrix.root(expectation$get()))
		},

		generate = function() {
			draw = matrix(rnorm(private$no.rows*nrow(private$expectation.root)), nrow=private$no.rows)
			if (!is.null(private$weights)) draw = sweep(draw, 1, private$weights, FUN="*")
			if (private$settings$exact) {
				if (nrow(draw) < ncol(draw)) input.error("cannot generate exact draws when number of rows is smaller than number of columns")
				draw = sweep(draw, 2, apply(draw, 2, mean), FUN="-")
				decomp = eigen(t(draw) %*% draw)
				draw = draw %*% decomp$vectors %*% diag.matrix(sqrt(1/decomp$values)) %*% private$expectation.root
				if (private$settings$mode == "covariance") draw = draw * sqrt(private$no.rows-1)
			} else {
				draw = draw %*% private$expectation.root
				if (private$settings$mode == "product") draw = draw / sqrt(private$no.rows)
			}
			return(draw)
		}
	)
)

