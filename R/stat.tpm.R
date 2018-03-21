#' Construct truncated product method statistic.
#' @param p - input p-values.
#' @param tau1 - truncation parameter. 0 < tau1 <= 1.
#' @return  Truncated product method statistic.
#' @details Let \eqn{p_{i}}, \eqn{i = 1,...,n} be a sequence of p-values, the TPM statistic
#' \deqn{TPM = \sum_{i=1}^n -2\log(p_i)I(p_i\leq\tau_2)}. TPM is the special case of TFisher when tau2=1.
#' @references 1. Hong Zhang and Zheyang Wu. "TFisher Tests: Optimal and Adaptive Thresholding for Combining p-Values", submitted.
#'
#' 2. Zaykin, D.V., Zhivotovsky, L. A., Westfall, P.H. and Weir, B.S. (2002), Truncated product method for combining P-values. Genet. Epidemiol., 22: 170–185. doi:10.1002/gepi.0042
#'
#' @examples
#' pval <- runif(100)
#' stat.tpm(p=pval, tau1=0.05)
#' @export

stat.tpm <- function(p,tau1){
  return(stat.tfisher(p,tau1,1))
}
