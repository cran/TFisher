#' Construct omnibus thresholding Fisher's p-value combination statistic.
#' @param p - input p-values.
#' @param TAU1 - a vector of truncation parameters. Must be in non-descending order.
#' @param TAU2 - a vector of normalization parameters. Must be in non-descending order.
#' @return  Omnibus thresholding Fisher's p-value combination statistic.
#' @details Let \eqn{p_{i}}, \eqn{i = 1,...,n} be a sequence of p-values, the thresholding statistics
#' \deqn{TFisher_j = \sum_{i=1}^n -2\log(p_i/\tau_{2j})I(p_i\leq\tau_{1j})}, \eqn{j = 1,...,d}.
#' The omnibus test statistic is the minimum p-value of these thresholding tests,
#' \deqn{W_o = min_j G_j(Soft_j)}, where \eqn{G_j} is the survival function of \eqn{Soft_j}.
#' @references 1. Hong Zhang and Zheyang Wu. "Optimal Thresholding of Fisher's P-value Combination
#' Tests for Signal Detection", submitted.
#'
#' @examples
#' pval = runif(20)
#' TAU1 = c(0.01, 0.05, 0.5, 1)
#' TAU2 = c(0.1, 0.2, 0.5, 1)
#' stat.tfisher.omni(p=pval, TAU1=TAU1, TAU2=TAU2)
#' @export

stat.tfisher.omni <- function(p, TAU1, TAU2){
  n = length(p)
  tfisherstat = mapply(function(x,y) stat.tfisher(p,x,y), TAU1, TAU2)
  pval = mapply(function(x,y,z) 1-p.tfisher(x, n, y, z), tfisherstat, TAU1, TAU2)
  return(min(pval))
}
