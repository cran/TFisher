#' Statistical power of thresholding Fisher's p-value combination test under Gaussian mixture model.
#' @param alpha - type-I error rate.
#' @param n - dimension parameter, i.e. the number of input p-values.
#' @param tau1 - truncation parameter. 0 < tau1 <= 1.
#' @param tau2 - normalization parameter. tau2 >= tau1.
#' @param eps - mixing parameter of the Gaussian mixture.
#' @param mu - mean of non standard Gaussian model.
#' @return Power of the thresholding Fisher's p-value combination test.
#' @details We consider the following hypothesis test,
#' \deqn{H_0: X_i\sim F_0, H_a: X_i\sim (1-\epsilon)F_0+\epsilon F_1}
#' , where \eqn{\epsilon} is the mixing parameter,
#' \eqn{F_0} is the standard normal CDF and \eqn{F = F_1} is the CDF of normal distribution with \eqn{\mu} defined by mu and \eqn{\sigma = 1}.
#'
#' @seealso \code{\link{stat.tfisher}} for the definition of the statistic.
#' @references 1. Hong Zhang and Zheyang Wu. "Optimal Thresholding of Fisher's P-value Combination
#' Tests for Signal Detection", submitted.
#' @examples
#' alpha = 0.05
#' #If the alternative hypothesis Gaussian mixture with eps = 0.1 and mu = 1.2:#
#' power.tfisher(alpha, 100, 0.05, 0.25, eps = 0.1, mu = 1.2)
#' @export
#' @importFrom sn psn
#' @importFrom stats pchisq integrate dbinom
#'
power.tfisher <- function(alpha,n,tau1,tau2,eps=0,mu=0){
  q = q.tfisher(p=1-alpha, n, tau1, tau2)
  if(eps==0|mu==0){
    return(1-(sum(pchisq(q+2*(1:n)*log(tau1/tau2),2*(1:n))*dbinom(1:n, n, tau1)) + (0<q)*dbinom(0, n, tau1)))
  }else{
    EY1 = integrate(function(x) -2*log(D_1(x, eps, mu)/tau2), 0, D(tau1,eps,mu))$value
    EY2 = integrate(function(x) (-2*log(D_1(x, eps, mu)/tau2))^2, 0, D(tau1,eps,mu))$value
    variance = n*(EY2-(EY1)^2)
    m = n*EY1
    EY3central = integrate(function(x) (-2*log(D_1(x, eps, mu)/tau2)-EY1)^3, 0, D(tau1,eps,mu), stop.on.error = FALSE)$value
    skewness = EY3central/(sqrt(n)*(EY2-(EY1)^2)^1.5)
    alpha = findsnalpha(skewness)
    omega = findsnomega(variance, alpha)
    xi = findsnxi(m, alpha, omega)
    return(1-psn(q,xi=xi,omega=omega,alpha=alpha))
  }
}
