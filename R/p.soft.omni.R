#' CDF of omnibus soft-thresholding Fisher's p-value combination statistic under the null hypothesis.
#' @param q - quantile, could be a vector.
#' @param n - dimension parameter, i.e. the number of p-values to be combined.
#' @param TAU1 - a vector of truncation parameters (=normalization parameters). Must be in non-descending order.
#' @return  The left-tail probability of the null distribution of omnibus soft-thresholding Fisher's p-value combination statistic.
#' @seealso \code{\link{stat.soft.omni}} for the definition of the statistic.
#' @references 1. Hong Zhang and Zheyang Wu. "Optimal Thresholding of Fisher's P-value Combination
#' Tests for Signal Detection", submitted.
#'
#' @examples
#' q = 0.01
#' n = 20
#' TAU1 = c(0.01, 0.05, 0.5, 1)
#' p.soft.omni(q=q, n=n, TAU1=TAU1)
#' @importFrom mvtnorm pmvnorm
#' @export

p.soft.omni <- function(q, n, TAU1){
  E = 2*n*TAU1
  nn = length(TAU1)
  M1 = matrix(rep(TAU1,nn),nn,nn)
  M2 = matrix(rep(2-TAU1+log(TAU1),each=nn,nn),nn,nn)
  COV = 4*n*(M1*M2-TAU1*log(TAU1))
  COV[lower.tri(COV)] = t(COV)[lower.tri(COV)]
  Q = mapply(function(x)q.soft(1-q, n, x),TAU1)

  #COV_half = t(solve(chol(COV)))
  #thr = COV_half%*%t(t(Q-E))
  #return(1-prod(pnorm(thr)))
  return(1-pmvnorm(lower=-Inf,upper=Q,mean=E,sigma=COV,abseps=1e-6)[1])
}
