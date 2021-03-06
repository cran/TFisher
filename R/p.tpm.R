#' CDF of truncated product method statistic under the null hypothesis.
#' @param q - quantile, could be a vector.
#' @param n - dimension parameter, i.e. the number of p-values to be combined.
#' @param tau1 - truncation parameter. 0 < tau1 <= 1.
#' @param M - correlation matrix of the input statistics. Default = NULL assumes independence.
#' @return The left-tail probability of the null distribution of truncated product method statistic at the given quantile.
#' @seealso \code{\link{stat.tpm}} for the definition of the statistic.
#' @references 1. Hong Zhang and Zheyang Wu. "TFisher Tests: Optimal and Adaptive Thresholding for Combining p-Values", submitted.
#'
#' 2. Zaykin, D.V., Zhivotovsky, L. A., Westfall, P.H. and Weir, B.S. (2002), Truncated product method for combining P-values. Genet. Epidemiol., 22: 170–185. doi:10.1002/gepi.0042
#'
#' @examples
#' pval <- runif(100)
#' tpmstat <- stat.tpm(p=pval, tau1=0.05)
#' p.tpm(q=tpmstat, n=100, tau1=0.05)
#' M = matrix(0.3,100,100) + diag(1-0.3,100)
#' p.tpm(q=tpmstat, n=100, tau1=0.05, M=M)
#' @export
#' @importFrom stats dbinom pchisq pgamma
#' @importFrom mvtnorm pmvnorm

p.tpm <- function(q,n,tau1,M=NULL){
  return(p.tfisher(q,n,tau1,1,M=M))
}
