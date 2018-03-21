#' Construct thresholding Fisher's p-value combination statistic.
#' @param p - input p-values.
#' @param tau1 - truncation parameter. 0 < tau1 <= 1.
#' @param tau2 - normalization parameter. tau2 >= tau1.
#' @return  Thresholding Fisher's p-value combination statistic.
#' @details Let \eqn{p_{i}}, \eqn{i = 1,...,n} be a sequence of p-values, the thresholding Fisher's p-value combination statistic
#' \deqn{TFisher = \sum_{i=1}^n -2\log(p_i/\tau_2)I(p_i\leq\tau_2)}
#' @references 1. Hong Zhang and Zheyang Wu. "TFisher Tests: Optimal and Adaptive Thresholding for Combining p-Values", submitted.
#' @examples
#' pval <- runif(100)
#' stat.tfisher(p=pval, tau1=0.05, tau2=0.25)
#' @export

stat.tfisher <- function(p,tau1,tau2){
  return(sum(-2*log(p[p<=tau1]/tau2)))
}

