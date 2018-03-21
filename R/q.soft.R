#' Quantile of soft-thresholding Fisher's p-value combination statistic under the null hypothesis.
#' @param p -  a scalar left probability that defines the quantile.
#' @param n - dimension parameter, i.e. the number of input p-values.
#' @param tau1 - truncation parameter=normalization parameter. tau1 > 0.
#' @param M - correlation matrix of the input statistics. Default = NULL assumes independence.
#' @return Quantile of soft-thresholding Fisher's p-value combination statistic.
#' @seealso \code{\link{stat.soft}} for the definition of the statistic.
#' @references 1. Hong Zhang and Zheyang Wu. "TFisher Tests: Optimal and Adaptive Thresholding for Combining p-Values", submitted.
#'
#' @examples
#' ## The 0.05 critical value of soft-thresholding statistic when n = 10:
#' q.soft(p=.99, n=20, tau1 = 0.05)
#' M = matrix(0.9,20,20) + diag(1-0.9,20)
#' q.soft(p=.99, n=20, tau1 = 0.05, M=M)
#' @export
#' @importFrom mvtnorm pmvnorm
#'
q.soft <- function(p, n, tau1,M=NULL){
  q.tfisher(p, n, tau1, tau1, M=M)
}
