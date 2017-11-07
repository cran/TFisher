#' CDF of omnibus thresholding Fisher's p-value combination statistic under the null hypothesis.
#' @param q - quantile, could be a vector.
#' @param n - dimension parameter, i.e. the number of p-values to be combined.
#' @param TAU1 - a vector of truncation parameters. Must be in non-descending order.
#' @param TAU2 - a vector of normalization parameters. Must be in non-descending order.
#' @return  The left-tail probability of the null distribution of omnibus thresholding Fisher's p-value combination statistic.
#' @seealso \code{\link{stat.tfisher.omni}} for the definition of the statistic.
#' @references 1. Hong Zhang and Zheyang Wu. "Optimal Thresholding of Fisher's P-value Combination
#' Tests for Signal Detection", submitted.
#'
#' @examples
#' q = 0.05
#' n = 20
#' TAU1 = c(0.01, 0.05, 0.5, 1)
#' TAU2 = c(0.1, 0.2, 0.5, 1)
#' p.tfisher.omni(q=q, n=n, TAU1=TAU1, TAU2=TAU2)
#' @importFrom mvtnorm pmvnorm
#' @export

p.tfisher.omni <- function(q, n, TAU1, TAU2){
  E = 2*n*TAU1*(1+log(TAU2/TAU1))
  nn = length(TAU1)

  minTAU = pmin(matrix(TAU1,nn,nn), t(matrix(TAU1,nn,nn)))
  part1 = 4*n*minTAU
  part2 = part1*(1+log(matrix(TAU2,nn,nn)/minTAU))*(1+log(t(matrix(TAU2,nn,nn))/minTAU))
  part3 = 4*n*matrix(TAU1,nn,nn)*t(matrix(TAU1,nn,nn))*(1+log(matrix(TAU2,nn,nn)/matrix(TAU1,nn,nn)))*(1+log(t(matrix(TAU2,nn,nn))/t(matrix(TAU1,nn,nn))))
  COV = part1 + part2 - part3

  COV[lower.tri(COV)] = t(COV)[lower.tri(COV)]
  Q = mapply(function(x,y)q.tfisher(1-q, n, x,y),TAU1,TAU2)

  #COV_half = t(solve(chol(COV)))
  #thr = COV_half%*%t(t(Q-E))
  #return(1-prod(pnorm(thr)))
  return(1-pmvnorm(lower=-Inf,upper=Q,mean=E,sigma=COV,abseps=1e-6)[1])
}
