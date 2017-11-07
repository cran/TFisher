#' Statistical power of truncated product method test under Gaussian mixture model.
#' @param alpha - type-I error rate.
#' @param n - dimension parameter, i.e. the number of input p-values.
#' @param tau1 - truncation parameter. 0 < tau1 <= 1. tau1 > 0.
#' @param eps - mixing parameter of the Gaussian mixture.
#' @param mu - mean of non standard Gaussian model.
#' @return Power of the truncated product method test.
#' @details We consider the following hypothesis test,
#' \deqn{H_0: X_i\sim F_0, H_a: X_i\sim (1-\epsilon)F_0+\epsilon F_1}
#' , where \eqn{\epsilon} is the mixing parameter,
#' \eqn{F_0} is the standard normal CDF and \eqn{F = F_1} is the CDF of normal distribution with \eqn{\mu} defined by mu and \eqn{\sigma = 1}.
#'
#' @seealso \code{\link{stat.soft}} for the definition of the statistic.
#' @references 1. Hong Zhang and Zheyang Wu. "Optimal Thresholding of Fisher's P-value Combination
#' Tests for Signal Detection", submitted.
#' @examples
#' alpha = 0.05
#' #If the alternative hypothesis Gaussian mixture with eps = 0.1 and mu = 1.2:#
#' power.tpm(alpha, 100, 0.05, eps = 0.1, mu = 1.2)
#' @export
#' @importFrom sn psn
power.tpm <- function(alpha,n,tau1,eps=0,mu=0){
  power.tfisher(alpha,n,tau1,1,eps,mu)
}
