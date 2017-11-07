#' CDF of thresholding Fisher's p-value combination statistic under the null hypothesis.
#' @param q - quantile, could be a vector.
#' @param n - dimension parameter, i.e. the number of p-values to be combined.
#' @param tau1 - truncation parameter. 0 < tau1 <= 1.
#' @param tau2 - normalization parameter. tau2 >= tau1.
#' @return The left-tail probability of the null distribution of thresholding Fisher's p-value combination statistic at the given quantile.
#' @seealso \code{\link{stat.tfisher}} for the definition of the statistic.
#' @references 1. Hong Zhang and Zheyang Wu. "Optimal Thresholding of Fisher's P-value Combination
#' Tests for Signal Detection", submitted.
#' @examples
#' pval <- runif(100)
#' tfstat <- stat.tfisher(p=pval, tau1=0.75, tau2=0.25)
#' p.tfisher(q=tfstat, n=100, tau1=0.75, tau2=0.25)
#' @export
#' @importFrom stats dbinom pchisq

p.tfisher <- function(q,n,tau1,tau2){
  return((sum(pchisq(q+2*(1:n)*log(tau1/tau2),2*(1:n))*dbinom(1:n, n, tau1)) + (0<q)*dbinom(0, n, tau1)))
}

