% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logistic.sim.R
\name{logistic.sim}
\alias{logistic.sim}
\title{Simulates time-series data from the Logistic map}
\usage{
logistic.sim(a = 4, s = 0, x0 = runif(1, 0, 1), n = 1000, n.start = 50)
}
\arguments{
\item{a}{a non-negative integer denoting the value of parameter \code{a} (Default 4).}

\item{s}{a non-negative integer denoting the variance value of the error term. If \eqn{s=0} gives the standard deterministic map (Default 0).}

\item{x0}{a non-negative integer denoting the initial condition (Default random number between 0 and 1).}

\item{n}{a non-negative integer denoting the length (Default 1000).}

\item{n.start}{a non-negative integer denoting the number of observations that will be discarded to ensure that the values are in the attractor (Default 50).}
}
\value{
A time-series data object generated from the Logistic map with or without an additive measurement noise term. This dataset could be useful for researchers interested in the field of chaotic dynamic systems and non-linear time series analysis and professors (and students) who teach (learn) courses related to those topics.
}
\description{
This function simulates time-series data from the Logistic map considering the parameter set selected by the user. The initial condition is a random number between 0 and 1. Some initial conditions may lead to an unstable system that will tend to infinity.
}
\note{
This function provides also noisy time-series data from the deterministic logistic map adding an additive measurement noise term if \eqn{s>0}. We have added to each time-series data a normal multinomial error term denoted by \eqn{{\varepsilon _t} \sim N\left( {0,s} \right)} with different variance values (\eqn{s}). In this sense we have considered it appropriate to add a measurement noise term because most real-world observed time-series data are usually noise-contaminated signals, characterised by an erratic and persistent volatility in certain periods and there is almost always a source of noise linked to measurement errors in real-world datasets.
}
\examples{
## set.seed(34)
## Simulates time-series data from the deterministic logistic map
## with a chaotic behaviour.
## ts <- logistic.sim(a=4, s=0, n=1000)
##
## Simulates time-series data from the deterministic logistic map
## with a non-chaotic behaviour.
## ts <- logistic.sim(a=3.2, s=0, n=1000)
}
\references{
May, R.M. 1976 Simple mathematical models with very complicated dynamics. Nature (261):459-467.
}
\author{
Julio E. Sandubete, Lorenzo Escot
}
