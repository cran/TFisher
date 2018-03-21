#' CDF of omnibus thresholding Fisher's p-value combination statistic under the null hypothesis.
#' @param q - quantile, could be a vector.
#' @param n - dimension parameter, i.e. the number of p-values to be combined.
#' @param TAU1 - a vector of truncation parameters. Must be in non-descending order.
#' @param TAU2 - a vector of normalization parameters. Must be in non-descending order.
#' @param M - correlation matrix of the input statistics. Default = NULL assumes independence.
#' @param P0 - a vector of point masses of TFisher statistics. Default = NULL.
#' @return  The left-tail probability of the null distribution of omnibus thresholding Fisher's p-value combination statistic.
#' @seealso \code{\link{stat.tfisher.omni}} for the definition of the statistic.
#' @references 1. Hong Zhang and Zheyang Wu. "TFisher Tests: Optimal and Adaptive Thresholding for Combining p-Values", submitted.
#' @examples
#' q = 0.05
#' n = 20
#' TAU1 = c(0.01, 0.05, 0.5, 1)
#' TAU2 = c(0.1, 0.2, 0.5, 1)
#' M = matrix(0.3,20,20) + diag(1-0.3,20)
#' p.tfisher.omni(q=q, n=n, TAU1=TAU1, TAU2=TAU2, M=M)
#' @export
#' @importFrom mvtnorm pmvnorm
#' @importFrom Matrix nearPD
#'
p.tfisher.omni <- function(q, n, TAU1, TAU2, M=NULL,P0=NULL){
  if(is.null(M)){
    E = 2*n*TAU1*(1+log(TAU2/TAU1))
    nn = length(TAU1)

    minTAU = pmin(matrix(TAU1,nn,nn), t(matrix(TAU1,nn,nn)))
    part1 = 4*n*minTAU
    part2 = part1*(1+log(matrix(TAU2,nn,nn)/minTAU))*(1+log(t(matrix(TAU2,nn,nn))/minTAU))
    part3 = 4*n*matrix(TAU1,nn,nn)*t(matrix(TAU1,nn,nn))*(1+log(matrix(TAU2,nn,nn)/matrix(TAU1,nn,nn)))*(1+log(t(matrix(TAU2,nn,nn))/t(matrix(TAU1,nn,nn))))
    COV = part1 + part2 - part3

    COV[lower.tri(COV)] = t(COV)[lower.tri(COV)]
    result = rep(NA, length(q))
    for(i in 1:length(q)){
      Q = mapply(function(x,y)q.tfisher(1-q, n, x,y),TAU1,TAU2)
      result[i] = 1-pmvnorm(lower=-Inf,upper=Q,mean=E,sigma=COV,abseps=1e-6)[1]
    }

    return(result)
  }else{
    p = length(TAU1)
    if(is.null(P0)){
      #truncation parameter
      bound = qnorm(1-TAU1/2)
      P0 = sapply(1:p,function(x)pmvnorm(lower=-rep(bound[x],n),upper=rep(bound[x],n),sigma=M)[1])
    }
    if(1-min(P0)<q){
      return("q is too large, impossible to obtain")
    }else{
      id = which(1-P0<q)
      if(length(id>0)){
        TAU1 = TAU1[-id]; TAU2 = TAU2[-id]
      }

      E = 2*n*TAU1*(1+log(TAU2/TAU1))

      COV = getTFisherCovM(M, TAU1, TAU2)
      if(any(eigen(COV)$value<0)){
        COV = as.matrix(nearPD(COV)$mat)
      }


      #Q = mapply(function(x,y)q.tfisher(1-q, n, x,y, M),TAU1,TAU2)
      Q = mapply(function(x,y)q.tfisher_N(1-q, n, x,y, M),TAU1,TAU2)
      result = 1-pmvnorm(lower=-Inf,upper=Q,mean=E,sigma=COV,abseps=1e-6)[1]
      return(result)
    }
  }
}
