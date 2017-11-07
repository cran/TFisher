#' Construct omnibus truncated product method statistic.
#' @param p - input p-values.
#' @param TAU1 - a vector of truncation parameters. Must be in non-descending order.
#' @return  Omnibus truncated product method statistic.
#' @details Let \eqn{p_{i}}, \eqn{i = 1,...,n} be a sequence of p-values, the truncated product method statistics
#' \deqn{TPM_j = \sum_{i=1}^n -2\log(p_i)I(p_i\leq\tau_{1j})}, \eqn{j = 1,...,d}.
#' The omnibus test statistic is the minimum p-value of these truncated product method tests,
#' \deqn{W_o = min_j G_j(TPM_j)}, where \eqn{G_j} is the survival function of \eqn{TPM_j}.
#' @references 1. Hong Zhang and Zheyang Wu. "Optimal Thresholding of Fisher's P-value Combination
#' Tests for Signal Detection", submitted.
#'
#' @examples
#' pval = runif(20)
#' TAU1 = c(0.01, 0.05, 0.5, 1)
#' stat.tpm.omni(p=pval, TAU1=TAU1)
#' @export

stat.tpm.omni <- function(p, TAU1){
  n = length(p)
  tpmstat = mapply(function(x) stat.tpm(p,x), TAU1)
  pval = mapply(function(x,y) 1-p.tpm(x, n, y), tpmstat, TAU1)
  return(min(pval))
}
