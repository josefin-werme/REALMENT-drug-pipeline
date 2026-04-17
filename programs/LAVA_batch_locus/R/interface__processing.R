

process.locus = function(data.set, locus.info, phenotypes=NULL, units=NULL, ..., catch.failures=T) {
	check.types(data.set="DataSet", locus.info="LocusDefinition")
	input = catchErrors(data.set$get.interface(phenotypes=phenotypes, ...)$set.locus(locus.info, detect.units=T)$set.units(units), all.failures=catch.failures)
	if (!is.null(input)) return(catchErrors(input$process(), all.failures=catch.failures))
}



#'
#'       TODO -> rewrite
#'
#'
#' Re-process locus to meta-analyse of selected phenotypes
#'
#' Will combine all elements of the requested phenotypes using standard inverse variance weighting, allowing them to be analysed as a single phenotype via the  multivariate analysis functions.
#' Note that the univariate test cannot currently be applied to meta-analysed phenotypes, so please do that beforehand on each phenotype individually.
#'
#' @param locus Locus object defined using the \code{\link{process.locus}} function.
#' @param meta.phenos Phenotypes you want to meta-analyse
#'
#' #'
#' @return This function returns an object just like that \code{\link{process.locus}} function, containing general locus info, the relevant processed sumstats, and info about the input phenotypes.
#'
#' \itemize{
#'     \item id - locus ID
#'     \item chr/start/stop - locus coordinates
#'     \item snps - list of locus SNPs
#'     \item N.snps - number of SNPs
#'     \item K - number of PCs
#'     \item delta - PC projected joint SNP effects for each phenotype
#'     \item sigma - sampling covariance matrix
#'     \item omega - genetic covariance matrix
#'     \item omega.cor - genetic correlation matrix
#'     \item N - vector of average N across locus SNPs for each phenotype
#'     \item phenos - phenotype IDs
#'     \item binary - boolean vector indicating whether phentoypes are binary
#' }
#'
#' @export


## for ProcessedLocus 'input' argument, create meta-analysed combination(s) of phenotypes (fixed effects, IVW) from vector(s) of phenotypes provided
## returns modified ProcessedLocus object with meta-analysed phenotypes added; component phenotypes used in meta-analysis are removed unless drop.components=F
## specification in ... should be in the form of [name]=[phenotypes], with [phenotypes] a vector of (partial) phenotype names to include; a meta-analyzed phenotype is created for each
meta.analyse = function(input, ..., drop.components=T) {
	check.types(input="ProcessedLocus")
	return(create.composite(input, ..., weights="ivw", drop.components=drop.components))
}



## for ProcessedLocus 'input' argument, create summed combination(s) of phenotypes from vector(s) of phenotypes provided
## returns modified ProcessedLocus object with composite phenotypes added; component phenotypes used for composite phenotypes are removed unless drop.components=F
## specification in ... should be in the form of [name]=[phenotypes], with [phenotypes] a vector of (partial) phenotype names to include; a meta-analyzed phenotype is created for each

## genetic component of composites are defined as weighted sums of input genetic components, scaled to a phenotypic variance of 1
## weights are relative, and are rescaled such that sum(abs(weights)) = 1; negative weights are permitted
## the weights argument can be specified as follows:
## - NULL or a single value: creates unweighted combination for all composite phenotypes
## - numeric vector: length must be equal to number of component phenotypes; can only be used when creating a single composite
## - list of numeric vectors: should contain an element for each composite phenotype (either a single value for unweighted, or a numeric vector of weights)
##   - order of weight specifications should be the same as the composite phenotype specifications, unless the list is named

## NB: composite phenotypes represent a combination of genetic components at the -genetic- level, and has no defined phenotypic variance
##     this will NOT correspond to the genetic component that would be obtained by taking the equivalent sum at the phenotype level and performing a GWAS, since the phenotypic correlations are unknown
##     consequently, only a h2.composite estimate is provided (and only if weights all have the same sign), which is the weighted mean of component h2 values
##     exception to this are meta-analysis composites, for which the phenotypic correlations are presumed to be one and phenotypic variance is therefore defined
create.composite = function(input, ..., weights=NULL, drop.components=F) {
	check.types(input="ProcessedLocus")
	phenotypes = list(...)
	if (length(phenotypes) == 0) input.error("no composite phenotype specification provided")
	if (is.null(names(phenotypes)) || any(names(phenotypes) == "")) input.error("missing names for composite phenotypes")
	IF.ANY(duplicate = duplicated(names(phenotypes)), THEN=input.error("duplicate composite phenotype %name% ", items.and=unique(names(phenotypes)[duplicate])))

	if (!is.list(weights)) {
		if (length(weights) > 1) weights = list(weights)
		else weights = as.list(rep(ifelse(length(weights) == 1, weights, 1), length(phenotypes)))
	} else if (!is.null(names(weights))) {
		IF.ANY(missing = !(names(phenotypes) %in% names(weights)), THEN=input.error("no weight %specification% provided for composite %phenotype%", items.and=names(phenotypes)[missing]))
		weights = weights[names(phenotypes)]
	}
	if (length(weights) != length(phenotypes)) input.error("number of weight specifications does not match number of composite phenotypes requested")

	return(input$add.composite(names=names(phenotypes), phenotypes=phenotypes, weights=weights, discard.input=drop.components))
}





