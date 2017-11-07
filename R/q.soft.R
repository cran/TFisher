#' Quantile of soft-thresholding Fisher's p-value combination statistic under the null hypothesis.
#' @param p -  a scalar left probability that defines the quantile.
#' @param n - dimension parameter, i.e. the number of input p-values.
#' @param tau1 - truncation parameter=normalization parameter. tau1 > 0.
#' @return Quantile of soft-thresholding Fisher's p-value combination statistic.
#' @seealso \code{\link{stat.soft}} for the definition of the statistic.
#' @references 1. Hong Zhang and Zheyang Wu. "Optimal Thresholding of Fisher's P-value Combination
#' Tests for Signal Detection", submitted.
#'
#' @examples
#' ## The 0.05 critical value of soft-thresholding statistic when n = 10:
#' q.soft(p=.99, n=10, tau1 = 0.05)
#' @export
q.soft <- function(p, n, tau1){
  q.tfisher(p, n, tau1, tau1)
}
