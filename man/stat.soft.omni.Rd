% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stat.soft.omni.R
\name{stat.soft.omni}
\alias{stat.soft.omni}
\title{Construct omnibus soft-thresholding Fisher's p-value combination statistic.}
\usage{
stat.soft.omni(p, TAU1, M = NULL)
}
\arguments{
\item{p}{- input p-values.}

\item{TAU1}{- a vector of truncation parameters (=normalization parameters). Must be in non-descending order.}

\item{M}{- correlation matrix of the input statistics. Default = NULL assumes independence.}
}
\value{
omni - omnibus soft-thresholding statistic.

pval - p-values of each soft-thresholding tests.
}
\description{
Construct omnibus soft-thresholding Fisher's p-value combination statistic.
}
\details{
Let \eqn{x_{i}}, \eqn{i = 1,...,n} be a sequence of individual statistics with
correlation matrix M, \eqn{p_{i}} be the corresponding two-sided p-values, then the soft-thresholding statistics
\deqn{Soft_j = \sum_{i=1}^n -2\log(p_i/\tau_{1j})I(p_i\leq\tau_{1j})}, \eqn{j = 1,...,d}.
The omnibus test statistic is the minimum p-value of these soft-thresholding tests,
\deqn{W_o = min_j G_j(Soft_j)}, where \eqn{G_j} is the survival function of \eqn{Soft_j}.
}
\examples{
pval = runif(20)
TAU1 = c(0.01, 0.05, 0.5, 1)
stat.soft.omni(p=pval, TAU1=TAU1)
M = matrix(0.3,20,20) + diag(1-0.3,20)
stat.soft.omni(p=pval, TAU1=TAU1, M=M)
}
\references{
1. Hong Zhang and Zheyang Wu. "TFisher Tests: Optimal and Adaptive Thresholding for Combining p-Values", submitted.
}
