% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select.R
\name{estProb}
\alias{estProb}
\title{Estimate the marginal posterior probability}
\usage{
estProb(iter_res, burnin)
}
\arguments{
\item{iter_res}{output from BMVS()}

\item{burnin}{burn-in iterations}
}
\value{
estimated marginal posterior probability
}
\description{
P(Z = 1)
}
\examples{
estProb(BMVS_res, burnin = 5000)
}
