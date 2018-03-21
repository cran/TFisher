#' CDF of soft-thresholding Fisher's p-value combination statistic under the null hypothesis.
#' @param q - quantile, could be a vector.
#' @param n - dimension parameter, i.e. the number of p-values to be combined.
#' @param tau1 - truncation parameter=normalization parameter. tau1 > 0.
#' @param M - correlation matrix of the input statistics. Default = NULL assumes independence.
#' @return The left-tail probability of the null distribution of soft-thresholding Fisher's p-value combination statistic at the given quantile.
#' @seealso \code{\link{stat.soft}} for the definition of the statistic.
#' @references 1. Hong Zhang and Zheyang Wu. "TFisher Tests: Optimal and Adaptive Thresholding for Combining p-Values", submitted.
#'
#' @examples
#' pval <- runif(100)
#' softstat <- stat.soft(p=pval, tau1=0.05)
#' p.soft(q=softstat, n=100, tau1=0.05)
#' M = matrix(0.3,100,100) + diag(1-0.3,100)
#' p.soft(q=softstat, n=100, tau1=0.05, M=M)
#' @export
#' @importFrom stats dbinom pchisq pgamma
#' @importFrom mvtnorm pmvnorm

p.soft <- function(q,n,tau1,M=NULL){
  return(p.tfisher(q,n,tau1,tau1,M=M))
}
