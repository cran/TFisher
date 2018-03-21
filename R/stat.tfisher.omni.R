#' Construct omnibus thresholding Fisher's (TFisher) p-value combination statistic.
#' @param p - input p-values from potentially correlated input sstatistics.
#' @param TAU1 - a vector of truncation parameters. Must be in non-descending order.
#' @param TAU2 - a vector of normalization parameters. Must be in non-descending order.
#' @param M - correlation matrix of the input statistics. Default = NULL assumes independence
#' @param MU - a vector of means of TFisher statistics. Default = NULL.
#' @param SIGMA2 - a vector of variances of TFisher statistics. Default = NULL.
#' @param P0 - a vector of point masses of TFisher statistics. Default = NULL.
#' @return  omni - omnibus TFisher statistic.
#' @return  pval - p-values of each TFisher tests.
#' @details Let \eqn{x_{i}}, \eqn{i = 1,...,n} be a sequence of individual statistics with
#' correlation matrix M, \eqn{p_{i}} be the corresponding two-sided p-values, then the TFisher statistics
#' \deqn{TFisher_j = \sum_{i=1}^n -2\log(p_i/\tau_{2j})I(p_i\leq\tau_{1j})}, \eqn{j = 1,...,d}.
#' The omnibus test statistic is the minimum p-value of these thresholding tests,
#' \deqn{W_o = min_j G_j(Soft_j)}, where \eqn{G_j} is the survival function of \eqn{Soft_j}.
#' @references 1. Hong Zhang and Zheyang Wu. "TFisher Tests: Optimal and Adaptive Thresholding for Combining p-Values", submitted.
#'
#' @examples
#' pval = runif(20)
#' TAU1 = c(0.01, 0.05, 0.5, 1)
#' TAU2 = c(0.1, 0.2, 0.5, 1)
#' stat.tfisher.omni(p=pval, TAU1=TAU1, TAU2=TAU2)
#' M = matrix(0.3,20,20) + diag(1-0.3,20)
#' stat.tfisher.omni(p=pval, TAU1=TAU1, TAU2=TAU2, M=M)
#' @export
#' @importFrom mvtnorm pmvnorm

stat.tfisher.omni <- function(p, TAU1, TAU2, M=NULL,MU=NULL,SIGMA2=NULL,P0=NULL){
  n = length(p)
  tfisherstat = mapply(function(x,y) stat.tfisher(p,x,y), TAU1, TAU2)
  if(!is.null(MU)&!is.null(SIGMA2)&!is.null(P0)){
    pvals = mapply(function(x,y,z,a,b,c) 1-p.tfisher(x, n, y, z, M,a,b,c), tfisherstat, TAU1, TAU2,MU,SIGMA2,P0)
  }else{
    pvals = mapply(function(x,y,z) 1-p.tfisher(x, n, y, z, M), tfisherstat, TAU1, TAU2)
  }
  return(list(omni=min(pvals),pvals=pvals))
}
