#' Construct soft-thresholding Fisher's p-value combination statistic.
#' @param p - input p-values.
#' @param tau1 - truncation parameter=normalization parameter. tau1 > 0.
#' @return  Soft-thresholding Fisher's p-value combination statistic.
#' @details Let \eqn{p_{i}}, \eqn{i = 1,...,n} be a sequence of p-values, the soft-thresholding statistic
#' \deqn{Soft = \sum_{i=1}^n -2\log(p_i/\tau_1)I(p_i\leq\tau_1)}. Soft-thresholding is the special case of TFisher when tau1=tau2.
#' @references 1. Hong Zhang and Zheyang Wu. "TFisher Tests: Optimal and Adaptive Thresholding for Combining p-Values", submitted.
#'
#' @examples
#' pval <- runif(100)
#' stat.soft(p=pval, tau1=0.05)
#' @export

stat.soft <- function(p,tau1){
  return(stat.tfisher(p,tau1,tau1))
}
