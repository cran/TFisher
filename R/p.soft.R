#' CDF of soft-thresholding Fisher's p-value combination statistic under the null hypothesis.
#' @param q - quantile, could be a vector.
#' @param n - dimension parameter, i.e. the number of p-values to be combined.
#' @param tau1 - truncation parameter=normalization parameter. tau1 > 0.
#' @return The left-tail probability of the null distribution of soft-thresholding Fisher's p-value combination statistic at the given quantile.
#' @seealso \code{\link{stat.soft}} for the definition of the statistic.
#' @references 1. Hong Zhang and Zheyang Wu. "Optimal Thresholding of Fisher's P-value Combination
#' Tests for Signal Detection", submitted.
#'
#' @examples
#' pval <- runif(100)
#' softstat <- stat.soft(p=pval, tau1=0.05)
#' p.soft(q=softstat, n=100, tau1=0.05)
#' @export

p.soft <- function(q,n,tau1){
  return(p.tfisher(q,n,tau1,tau1))
}
