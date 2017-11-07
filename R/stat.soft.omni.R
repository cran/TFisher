#' Construct omnibus soft-thresholding Fisher's p-value combination statistic.
#' @param p - input p-values.
#' @param TAU1 - a vector of truncation parameters (=normalization parameters). Must be in non-descending order.
#' @return  Omnibus soft-thresholding Fisher's p-value combination statistic.
#' @details Let \eqn{p_{i}}, \eqn{i = 1,...,n} be a sequence of p-values, the soft-thresholding statistics
#' \deqn{Soft_j = \sum_{i=1}^n -2\log(p_i/\tau_{1j})I(p_i\leq\tau_{1j})}, \eqn{j = 1,...,d}.
#' The omnibus test statistic is the minimum p-value of these soft-thresholding tests,
#' \deqn{W_o = min_j G_j(Soft_j)}, where \eqn{G_j} is the survival function of \eqn{Soft_j}.
#' @references 1. Hong Zhang and Zheyang Wu. "Optimal Thresholding of Fisher's P-value Combination
#' Tests for Signal Detection", submitted.
#'
#' @examples
#' pval = runif(20)
#' TAU1 = c(0.01, 0.05, 0.5, 1)
#' stat.soft.omni(p=pval, TAU1=TAU1)
#' @export

stat.soft.omni <- function(p, TAU1){
  n = length(p)
  softstat = mapply(function(x) stat.soft(p,x), TAU1)
  pval = mapply(function(x,y) 1-p.soft(x, n, y), softstat, TAU1)
  return(min(pval))
}
