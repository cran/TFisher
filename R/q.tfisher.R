#' Quantile of thresholding Fisher's p-value combination statistic under the null hypothesis.
#' @param p -  a scalar left probability that defines the quantile.
#' @param n - dimension parameter, i.e. the number of input p-values.
#' @param tau1 - truncation parameter. 0 < tau1 <= 1.
#' @param tau2 - normalization parameter. tau2 >= tau1.
#' @return Quantile of thresholding Fisher's p-value combination statistic.
#' @seealso \code{\link{stat.tfisher}} for the definition of the statistic.
#' @references 1. Hong Zhang and Zheyang Wu. "Optimal Thresholding of Fisher's P-value Combination
#' Tests for Signal Detection", submitted.
#' @examples
#' ## The 0.05 critical value of TFisher statistic when n = 10:
#' q.tfisher(p=.95, n=10, tau1 = 0.05, tau2 = 0.25)
#' @export
#'
q.tfisher <- function(p, n, tau1, tau2){
  if(p<dbinom(0, n, tau1)){
    return(0)
  }
  lower = 0
  upper = max(2*n*(1+log(tau2/tau1))+6*sqrt(4*n*tau1*(1+(1-tau1)*(1+log(tau2/tau1))^2)))
  maxnumIter = 500
  numIter = 1
  p_cal_lower = p.tfisher(lower, n, tau1, tau2)
  p_cal_upper = p.tfisher(upper, n, tau1, tau2)

  while ((p_cal_lower - p) * (p_cal_upper - p) > 0 && numIter <
         maxnumIter) {
    if (p_cal_lower > p) {
      lower = lower * 1.1
      p_cal_lower = p.tfisher(lower, n, tau1, tau2)
    }
    if (p_cal_upper < p) {
      upper = upper * 1.1
      p_cal_upper = p.tfisher(upper, n, tau1, tau2)
    }
    numIter = numIter + 1
  }
  q = mean(c(lower, upper))
  p_cal = p.tfisher(q, n, tau1, tau2)
  error = (p_cal - p)/p
  numIter = 1
  while (abs(error) > 1e-08 && numIter < maxnumIter) {
    if (error > 0) {
      upper = q
      q = mean(c(lower, upper))
      p_cal = p.tfisher(q, n, tau1, tau2)
      error = (p_cal - p)/p
    }
    else {
      lower = q
      q = mean(c(lower, upper))
      p_cal = p.tfisher(q, n, tau1, tau2)
      error = (p_cal - p)/p
    }
    numIter = numIter + 1
  }
  if (numIter < maxnumIter) {
    return(q)
  }
  else {
    return(paste0("The quantile is larger than ", upper))
  }
}
