#' CDF of omnibus truncated product method statistic under the null hypothesis.
#' @param q - quantile, could be a vector.
#' @param n - dimension parameter, i.e. the number of p-values to be combined.
#' @param TAU1 - a vector of truncation parameters. Must be in non-descending order.
#' @return  The left-tail probability of the null distribution of omnibus truncated product method statistic.
#' @seealso \code{\link{stat.tpm.omni}} for the definition of the statistic.
#' @references 1. Hong Zhang and Zheyang Wu. "Optimal Thresholding of Fisher's P-value Combination
#' Tests for Signal Detection", submitted.
#'
#' @examples
#' q = 0.05
#' n = 20
#' TAU1 = c(0.01, 0.05, 0.5, 1)
#' p.tpm.omni(q=q, n=n, TAU1=TAU1)
#' @importFrom mvtnorm pmvnorm
#' @export

p.tpm.omni <- function(q, n, TAU1){
  p.tfisher.omni(q, n, TAU1, rep(1,length(TAU1)))
}
