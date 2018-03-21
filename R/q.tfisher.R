#' Quantile of thresholding Fisher's p-value combination statistic under the null hypothesis.
#' @param p -  a scalar left probability that defines the quantile.
#' @param n - dimension parameter, i.e. the number of input p-values.
#' @param tau1 - truncation parameter. 0 < tau1 <= 1.
#' @param tau2 - normalization parameter. tau2 >= tau1.
#' @param M - correlation matrix of the input statistics. Default = NULL assumes independence.
#' @return Quantile of thresholding Fisher's p-value combination statistic.
#' @seealso \code{\link{stat.tfisher}} for the definition of the statistic.
#' @references 1. Hong Zhang and Zheyang Wu. "TFisher Tests: Optimal and Adaptive Thresholding for Combining p-Values", submitted.
#' @examples
#' ## The 0.05 critical value of TFisher statistic when n = 10:
#' q.tfisher(p=.95, n=20, tau1=0.05, tau2=0.25)
#' ## when corrrelated
#' M = matrix(0.3,20,20) + diag(1-0.3,20)
#' q.tfisher(p=.95, n=20, tau1=0.05, tau2=0.25, M=M)
#' @export
#' @importFrom mvtnorm pmvnorm
#'
q.tfisher <- function(p, n, tau1, tau2, M=NULL){
  if(is.null(M)){
    if(p<dbinom(0, n, tau1)){
      return(0)
    }else{
      lower = 0
      upper = max(2*n*(1+log(tau2/tau1))+6*sqrt(4*n*tau1*(1+(1-tau1)*(1+log(tau2/tau1))^2)))
      maxnumIter = 500
      numIter = 1
      p_cal_lower = p.tfisher(lower, n, tau1, tau2)
      p_cal_upper = p.tfisher(upper, n, tau1, tau2)

      while((p_cal_lower - p) * (p_cal_upper - p) > 0 && numIter <
            maxnumIter) {
        if(p_cal_lower > p) {
          lower = lower * 1.1
          p_cal_lower = p.tfisher(lower, n, tau1, tau2)
        }
        if(p_cal_upper < p) {
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
        if(error > 0) {
          upper = q
          q = mean(c(lower, upper))
          p_cal = p.tfisher(q, n, tau1, tau2)
          error = (p_cal - p)/p
        }else{
          lower = q
          q = mean(c(lower, upper))
          p_cal = p.tfisher(q, n, tau1, tau2)
          error = (p_cal - p)/p
        }
        numIter = numIter + 1
      }
      if(numIter < maxnumIter) {
        return(q)
      }else{
        return(paste0("The quantile is larger than ", upper))
      }
    }
  }else{
    n = length(M[1,])
    #truncation parameter
    bound = qnorm(1-tau1/2)
    p0 = pmvnorm(lower=-rep(bound,n),upper=rep(bound,n),sigma=M)[1]
    if(p<p0){
      return(0)
    }else{
      mu = getTFisherMean(n, tau1, tau2)
      sigma2 = getTFisherVar_v1(M, tau1, tau2)

      #gap parameter
      shift = -2*log(tau1/tau2)

      #gamma parameters
      a = (mu-shift*(1-p0))^2/((1-p0)*sigma2-p0*mu^2)
      s = (mu/(1-p0)-shift)/a
      Fw = function(q)(0<q)*p0+(1-p0)*pgamma(q-shift, shape=a, scale=s)
      lower = 0
      upper = max(2*n*(1+log(tau2/tau1))+6*sqrt(4*n*tau1*(1+(1-tau1)*(1+log(tau2/tau1))^2)))
      maxnumIter = 500
      numIter = 1
      p_cal_lower = Fw(lower)
      p_cal_upper = Fw(upper)

      while((p_cal_lower - p) * (p_cal_upper - p) > 0 && numIter <
            maxnumIter) {
        if(p_cal_lower > p) {
          lower = lower * 1.1
          p_cal_lower = Fw(lower)
        }
        if(p_cal_upper < p) {
          upper = upper * 1.1
          p_cal_upper = Fw(upper)
        }
        numIter = numIter + 1
      }
      q = mean(c(lower, upper))
      p_cal = Fw(q)
      error = (p_cal - p)/p
      numIter = 1
      while (abs(error) > 1e-08 && numIter < maxnumIter) {
        if(error > 0) {
          upper = q
          q = mean(c(lower, upper))
          p_cal = Fw(q)
          error = (p_cal - p)/p
        }else{
          lower = q
          q = mean(c(lower, upper))
          p_cal = Fw(q)
          error = (p_cal - p)/p
        }
        numIter = numIter + 1
      }
      if(numIter < maxnumIter) {
        return(q)
      }else{
        return(paste0("The quantile is larger than ", upper))
      }
    }
  }
}
