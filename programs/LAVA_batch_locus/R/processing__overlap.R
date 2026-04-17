#' @include global.R
#' @include utils__misc.R
#' @include processing__ld_blocks.R
#' @include processing__univariate.R



SamplingCorrelation = R6::R6Class("SamplingCorrelation",
	inherit = CovarianceMatrix,
	private = list(

	),
	public = list(

	)
)



compute.scorr = function(marg1, marg2, p.thresh=NULL, p.id=1:2) {
	# check same LD block, not truncated etc.

	delta = cbind(marg1$get("delta"), marg2$get("delta"))
	use = rep(T, nrow(delta))
	if (!is.null(p.thresh)) {
		sigma = c(marg1$get("sigma"), marg2$get("sigma"))
		pval = pchisq(sweep(delta^2, 2, sigma, FUN="/"), df=1, lower.tail=F)
		for (p in p.id) use = use & pval[,p] > p.thresh
	}


	ld = marg1$get.ld()
	Q = ld$get.vectors()
	L = ld$get.values()
	M = sqrt(L) * as.numeric(matrix(1, nrow=1, ncol=ld$no.snps()) %*% Q)

	x = cbind(1, L, M^2)[use,]
	y = cbind(delta^2, delta[,1]*delta[,2])[use,]

	corr = c()
	for (p in 1:3) {
		curr = x[,1:p,drop=F]
		est = (solve(t(curr) %*% curr) %*% t(curr) %*% y)[1,]
		corr[p] = est[3] / sqrt(prod(est[1:2])) # add negative check
	}
	return(corr)
}

