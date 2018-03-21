#' CDF of omnibus truncated product method statistic under the null hypothesis.
#' @param q - quantile, could be a vector.
#' @param n - dimension parameter, i.e. the number of p-values to be combined.
#' @param TAU1 - a vector of truncation parameters. Must be in non-descending order.
#' @param M - correlation matrix of the input statistics. Default = NULL assumes independence
#' @return  The left-tail probability of the null distribution of omnibus truncated product method statistic.
#' @seealso \code{\link{stat.tpm.omni}} for the definition of the statistic.
#' @references 1. Hong Zhang and Zheyang Wu. "TFisher Tests: Optimal and Adaptive Thresholding for Combining p-Values", submitted.
#'
#' @examples
#' q = 0.05
#' n = 20
#' TAU1 = c(0.01, 0.05, 0.5, 1)
#' M = matrix(0.3,20,20) + diag(1-0.3,20)
#' p.tpm.omni(q=q, n=n, TAU1=TAU1, M=M)
#' @export
#' @importFrom mvtnorm pmvnorm
#'
p.tpm.omni <- function(q, n, TAU1, M=NULL){
  p.tfisher.omni(q, n, TAU1, rep(1,length(TAU1)), M=M)
}
