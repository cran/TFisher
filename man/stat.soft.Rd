% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stat.soft.R
\name{stat.soft}
\alias{stat.soft}
\title{Construct soft-thresholding Fisher's p-value combination statistic.}
\usage{
stat.soft(p, tau1)
}
\arguments{
\item{p}{- input p-values.}

\item{tau1}{- truncation parameter=normalization parameter. tau1 > 0.}
}
\value{
Soft-thresholding Fisher's p-value combination statistic.
}
\description{
Construct soft-thresholding Fisher's p-value combination statistic.
}
\details{
Let \eqn{p_{i}}, \eqn{i = 1,...,n} be a sequence of p-values, the soft-thresholding statistic
\deqn{Soft = \sum_{i=1}^n -2\log(p_i/\tau_1)I(p_i\leq\tau_1)}. Soft-thresholding is the special case of TFisher when tau1=tau2.
}
\examples{
pval <- runif(100)
stat.soft(p=pval, tau1=0.05)
}
\references{
1. Hong Zhang and Zheyang Wu. "TFisher Tests: Optimal and Adaptive Thresholding for Combining p-Values", submitted.
}
