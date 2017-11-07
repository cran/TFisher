#' Quantile of truncated product method statistic under the null hypothesis.
#' @param p -  a scalar left probability that defines the quantile.
#' @param n - dimension parameter, i.e. the number of input p-values.
#' @param tau1 - truncation parameter. 0 < tau1 <= 1.
#' @return Quantile of truncated product method statistic.
#' @seealso \code{\link{stat.tpm}} for the definition of the statistic.
#' @references 1. Hong Zhang and Zheyang Wu. "Optimal Thresholding of Fisher's P-value Combination
#' Tests for Signal Detection", submitted.
#'
#' 2. Zaykin, D.V., Zhivotovsky, L. A., Westfall, P.H. and Weir, B.S. (2002), Truncated product method for combining P-values. Genet. Epidemiol., 22: 170–185. doi:10.1002/gepi.0042
#'
#' @examples
#' ## The 0.05 critical value of TPM statistic when n = 10:
#' q.tpm(p=.95, n=10, tau1 = 0.05)
#' @export
q.tpm <- function(p, n, tau1){
  q.tfisher(p, n, tau1, 1)
}
