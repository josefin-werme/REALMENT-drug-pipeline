

run.univ = function(locus.data, phenotypes=NULL) {
	check.types(locus.data="ProcessedLocus")
	if (!is.null(phenotypes)) locus.data = locus.data$subset(phenotypes, empty.mode="quiet")
	return(locus.data$summarize())
}


#'
#'       TODO -> rewrite
#'
#'
#' Bivariate local genetic correlation analysis
#'
#' Performs bivariate local genetic correlation analysis between two phenotypes.
#' By default, the bivariate test will be performed for all combinations of phenotypes in the locus,
#' but this can be modified using the 'phenos' and 'target' arguments (see below)
#'
#' @param locus Locus object created using the the \code{\link{process.locus}} function. Contains all the relevant parameters and processed sum-stats for the phenotypes of interest
#' @param phenos Subset of phenotypes to analyse. If NULL, all phenotypes in the locus object will be analysed
#' @param target Target phenotype of interest. If NULL, bivariate correlations between all pairs of phenotypes will be computed;
#' Otherwise, only the relations between the target phenotype and the other phenotypes will be tested.
#' @param adap.thresh The thresholds at which to increase the number of iterations for the p-value generation.
#' Default number of iterations is 1e+4, but will be increased to 1e+5, and 1e+6 as p-values fall below the respective thresholds.
#' If set to NULL, the maximum number of iterations is capped at the default (Note: this significantly speeds up the analysis, but results in poor accuracy for low p-values)
#' @param p.values Set to F to suppress p-values
#' @param CIs Set to F to suppress 95\% confidence intervals
#' @param param.lim The +- threshold at which estimated parameters are considered to be too far out of bounds. If the estimated parameter exceeds this threshold, it is considered unreliable and will be set to NA.
#'
#' @return Data frame with the columns:
#' \itemize{
#'     \item phen1 / phen2 - analysed phenotypes
#'     \item rho - standardised coefficient for the local genetic correlation
#'     \item rho.lower / rho.upper - 95\% confidence intervals for rho
#'     \item r2 - proportion of variance in genetic signal for phen1 explained by phen2 (and vice versa)
#'     \item r2.lower / r2.upper - 95\% confidence intervals for the r2
#'     \item p - simulation p-values for the local genetic correlation
#' }
#' @export



## NB: setting ci.threshold explicitly will automatically switch compute.ci="significant"


run.bivar = function(locus.data, targets=NULL, phenotypes=NULL, include.univariate=F, ..., catch.failures=T) {
	check.types(locus.data="ProcessedLocus")

	results = catchErrors(AnalysisProcessor$new(locus.data, ...)$bivariate(targets=targets, phenotypes=phenotypes), all.failures=catch.failures)
	if (include.univariate && !is.null(results) ) {
		phenotypes = unique(unlist(results[c("phenotype1", "phenotype2")]))
		return(list(univariate = locus.data$subset(phenotypes, empty.mode="quiet")$summarize(),	bivariate = results))
	} else return(results)
}



#'
#'       TODO -> rewrite
#'
#'
#' Local genetic multiple regression analysis
#'
#' Will perform a local genetic multiple regression analysis, which models the genetic signal for a single outcome phenotype of interest using two or more predictor phenotypes.
#' Here, the genetic correlations between all predictors will be accounted for, and their genetic relation with the outcome will be conditioned on one another.
#'
#' @param locus Locus object created using the the \code{\link{process.locus}} function. Contains all the relevant parameters and processed sum-stats for the phenotypes of interest
#' @param target Outcome phenotype of interest (all other phenotypes will be considered predictors)
#' @param phenos Subset of phenotypes to analyse. If NULL, all phenotypes in the locus object will be analysed.
#' @param adap.thresh The thresholds at which to increase the number of iterations for the p-value generation.
#' Default number of iterations is 1e+4, but will be increased to 1e+5, and 1e+6 as p-values fall below the respective thresholds.
#' If set to NULL, the maximum number of iterations is capped at the default (Note: this significantly speeds up the analysis, but results in poor accuracy for low p-values)
#' @param p.values Set to F to suppress p-values
#' @param CIs Set to F to suppress 95\% confidence intervals
#' @param param.lim The +- threshold at which estimated parameters are considered to be too far out of bounds. If the estimated parameter exceeds this threshold, it is considered unreliable and will be set to NA.
#'
#' @return Data frame with the columns:
#' \itemize{
#'     \item predictors / outcome - analysed phenotypes
#'     \item gamma - standardised multiple regression coefficient
#'     \item gamma.lower / gamma.upper - 95\% confidence intervals for gamma
#'     \item r2 - proportion of variance in genetic signal for the outcome phenotype explained by all predictor phenotypes simultaneously
#'     \item r2.lower / r2.upper - 95\% confidence intervals for the r2
#'     \item p - simulation p-values for the gammas
#' }
#' @export
#'
#'
#'



run.multireg = function(locus.data, outcomes, predictors=NULL, include.marginal=F, include.submodels=F, ..., catch.failures=T) {
	check.types(locus.data="ProcessedLocus")

	processor = AnalysisProcessor$new(locus.data, ...)
	results = catchErrors(processor$regression(outcomes=outcomes, predictors=predictors), all.failures=catch.failures)
	if (!is.null(results)) {
		if (include.marginal) {
			phenotypes = unique(c(results$models$outcome, unlist(results$models$predictors)))
			results$marginal = list(
				univariate = processor$get.univariate(phenotypes),
				rg = processor$get.rg(phenotypes)
			)
		}

		if (include.submodels && any(results$models$no.predictors > 1)) {
			sub.results = list()
			for (i in which(results$models$no.predictors > 1)) {
				model.name = results$models$model.name[i]
				sub.results = list.append(sub.results, list(models = results$models[i,], coefficients = results$coefficients[results$coefficients$model.name == model.name,]))
				subsets = make.subsets(results$models$predictors[[i]]); subsets = subsets[sapply(subsets, length) < results$models$no.predictors[i]]
				for (sub.pred in subsets) {
					curr = catchErrors(processor$regression(outcomes=results$models$outcome[i], predictors=sub.pred), all.failures=T, all.errors=T)
					if (!is.null(curr)) sub.results = list.append(sub.results, curr)
				}
			}
			sub.results = lapply(list.extract(sub.results, "models", "coefficients"), fill.rowmerge)
			sub.results$models = sub.results$models[order(sub.results$models$outcome, -sub.results$models$no.predictors, sub.results$models$model.name),]
			sub.results$coefficients = sub.results$coefficients[order(match(sub.results$coefficients$model.name, sub.results$models$model.name), sub.results$coefficients$predictor),]
			results$sub.models = sub.results
		}
	}
	return(results)
}





#'
#'       TODO -> rewrite
#'
#'
#' Local partial genetic correlation analysis
#'
#' Will perform a local partial genetic correlation between the first two phenotypes (phen1, phen2) conditioned on the rest (Z).
#' Phenotype order is based on that within the locus object by default, but can be changed by passing a phenotype vector with the desired order to the 'phenos' argument.
#'
#' @param locus Locus object created using the the \code{\link{process.locus}} function. Contains all the relevant parameters and processed sum-stats for the phenotypes of interest
#' @param target The two target phenotypes of interest for which the partial correlation will be computed. All other phenotypes will be conditioned on.
#' @param phenos Subset of phenotypes to analyse. If NULL, all phenotypes in the locus object will be analysed.
#' @param adapt.thresh The thresholds at which to increase the number of iterations for the p-value generation.
#' Default number of iterations is 1e+4, but will be increased to 1e+5, and 1e+6 as p-values fall below the respective thresholds.
#' If set to NULL, the maximum number of iterations is capped at the default (Note: this significantly speeds up the analysis, but results in poor accuracy for low p-values)
#' @param p.values Set to F to suppress p-values
#' @param CIs Set to F to suppress 95\% confidence intervals
#' @param max.r2 Max r2 threshold for the regression of phen1~Z and phen2~Z. If any of these r2's are too high, the partial correlation becomes unstable, and analysis is therefore aborted.
#' @param param.lim The +- threshold at which estimated parameters are considered to be too far out of bounds. If the estimated parameter exceeds this threshold, it is considered unreliable and will be set to NA.
#'
#' @return Data frame with the columns:
#' \itemize{
#'     \item phen1 / phen2 - target phenotypes
#'     \item z - phenotype(s) which the genetic correlation between the target phenotypes were conditioned on
#'     \item r2.phen1_z / r2.phen2_z - the proportion of genetic signal in target phenotypes explained by z. Note: if either of these exceed the threshold specified by the max.r2 argument, the analysis will be aborted (to prevent unstable estimates)
#'     \item pcor - the partial correlation between phen1 and phen2 conditioned on z
#'     \item ci.lower / ci.upper - 95\% confidence intervals for the partial genetic correlation
#'     \item p - simulation p-values for the partial genetic correlation
#' }
#' @export




run.pcor = function(locus.data, outcomes=NULL, predictors=NULL, include.marginal=F, include.submodels=F, ..., catch.failures=T) {
	check.types(locus.data="ProcessedLocus")

	processor = AnalysisProcessor$new(locus.data, ...)
	results = catchErrors(processor$partial(outcomes=outcomes, predictors=predictors), all.failures=catch.failures)
	if (!is.null(results)) {
		include = list()
		if (include.marginal) {
			phenotypes = unique(unlist(results[c("outcome1", "outcome2", "predictors")]))
			include$marginal = list(
				univariate = processor$get.univariate(phenotypes),
				rg = processor$get.rg(phenotypes)
			)
		}

		if (include.submodels && any(results$no.predictors > 1)) {
			sub.results = list()
			for (i in which(results$no.predictors > 1)) {
				sub.results = list.append(sub.results, results[i,])
				subsets = make.subsets(results$predictors[[i]]); subsets = subsets[sapply(subsets, length) < results$no.predictors[i]]
				for (sub.pred in subsets) {
					curr = catchErrors(processor$partial(outcomes=results[i,c("outcome1", "outcome2")], predictors=sub.pred), all.failures=T, all.errors=T)
					if (!is.null(curr)) sub.results = list.append(sub.results, curr)
				}
			}
			sub.results = fill.rowmerge(sub.results)
			include$sub.models = sub.results[order(sub.results$outcome1, sub.results$outcome2, -sub.results$no.predictors, sub.results$model.name),]
		}
		if (length(include) > 0) results = c(list(partial.rg = results), include)
	}
	return(results)
}









