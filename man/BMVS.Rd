% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BMVS.R
\name{BMVS}
\alias{BMVS}
\title{BMVS main function}
\usage{
BMVS(Y, X, iteration = 3000)
}
\arguments{
\item{Y}{input response variables}

\item{X}{input predictors}

\item{iteration}{the number of iterations of the MCMC chain}
}
\value{
The MCMC simulation resutls. The results can be analyzed by estProb(), estBeta() and select().
}
\description{
This is the main function of BMVS algorithm
}
\examples{
data(tobacco)
Y = as.matrix(tobacco[, 1:3])
X = as.matrix(tobacco[, 4:9])
BMVS_res = BMVS(Y, X, iteration = 10000)
}
