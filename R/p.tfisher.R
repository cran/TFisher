#' CDF of thresholding Fisher's p-value combination statistic under the null hypothesis.
#' @param q - quantile, could be a vector.
#' @param n - dimension parameter, i.e. the number of p-values to be combined.
#' @param tau1 - truncation parameter. 0 < tau1 <= 1.
#' @param tau2 - normalization parameter. tau2 >= tau1.
#' @param M - correlation matrix of the input statistics. Default = NULL assumes independence.
#' @param mu - the mean of TFisher statistics. Default = NULL.
#' @param sigma2 - the variance of TFisher statistics. Default = NULL.
#' @param p0 - the point masse of TFisher statistics. Default = NULL.
#' @return The left-tail probability of the null distribution of thresholding Fisher's p-value combination statistic at the given quantile.
#' @seealso \code{\link{stat.tfisher}} for the definition of the statistic.
#' @references 1. Hong Zhang and Zheyang Wu. "TFisher Tests: Optimal and Adaptive Thresholding for Combining p-Values", submitted.
#' @examples
#' pval <- runif(20)
#' tfstat <- stat.tfisher(p=pval, tau1=0.25, tau2=0.75)
#' p.tfisher(q=tfstat, n=20, tau1=0.25, tau2=0.75)
#' M = matrix(0.3,20,20) + diag(1-0.3,20)
#' p.tfisher(q=tfstat, n=20, tau1=0.25, tau2=0.75, M=M)
#' @export
#' @importFrom stats dbinom pchisq pgamma
#' @importFrom mvtnorm pmvnorm

p.tfisher <- function(q,n,tau1,tau2,M=NULL,mu=NULL,sigma2=NULL, p0=NULL){
  if(is.null(M)){
    return((sum(pchisq(q+2*(1:n)*log(tau1/tau2),2*(1:n))*dbinom(1:n, n, tau1)) + (0<q)*dbinom(0, n, tau1)))
  }else{
    n = length(M[1,])
    if(is.null(mu)){
      mu = getTFisherMean(n, tau1, tau2)
    }
    if(is.null(sigma2)){
      sigma2 = getTFisherVar_v1(M, tau1, tau2)
    }
    if(is.null(p0)){
      #truncation parameter
      bound = qnorm(1-tau1/2)
      p0 = pmvnorm(lower=-rep(bound,n),upper=rep(bound,n),sigma=M)[1]
    }


    #gap parameter
    shift = -2*log(tau1/tau2)

    #gamma parameters
    a = (mu-shift*(1-p0))^2/((1-p0)*sigma2-p0*mu^2)
    s = (mu/(1-p0)-shift)/a

    return((0<q)*p0+(1-p0)*pgamma(q-shift, shape=a, scale=s))
  }
}


