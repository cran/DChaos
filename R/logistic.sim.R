################################################################################
############ DChaos: An R Package for Chaotic Time Series Analysis #############
################################################################################
# FUNCTION:         DESCRIPTION:
#   logistic.sim      Simulates time series from the Logistic map
################################################################################
#' Simulates time-series data from the Logistic map
#' @name logistic.sim
#' @aliases logistic.sim
#' @description
#' This function simulates time-series data from the Logistic map considering the parameter set selected by the user. The initial condition is a random number between 0 and 1. Some initial conditions may lead to an unstable system that will tend to infinity.
#' @param a a non-negative integer denoting the value of parameter \code{a} (Default 4).
#' @param s a non-negative integer denoting the variance value of the error term. If \eqn{s=0} gives the standard deterministic map (Default 0).
#' @param x0 a non-negative integer denoting the initial condition (Default random number between 0 and 1).
#' @param n a non-negative integer denoting the length (Default 1000).
#' @param n.start a non-negative integer denoting the number of observations that will be discarded to ensure that the values are in the attractor (Default 50).
#' @return A time-series data object generated from the Logistic map with or without an additive measurement noise term. This dataset could be useful for researchers interested in the field of chaotic dynamic systems and non-linear time series analysis and professors (and students) who teach (learn) courses related to those topics.
#' @note This function provides also noisy time-series data from the deterministic logistic map adding an additive measurement noise term if \eqn{s>0}. We have added to each time-series data a normal multinomial error term denoted by \eqn{{\varepsilon _t} \sim N\left( {0,s} \right)} with different variance values (\eqn{s}). In this sense we have considered it appropriate to add a measurement noise term because most real-world observed time-series data are usually noise-contaminated signals, characterised by an erratic and persistent volatility in certain periods and there is almost always a source of noise linked to measurement errors in real-world datasets.
#' @references May, R.M. 1976 Simple mathematical models with very complicated dynamics. Nature (261):459-467.
#' @examples
#' ## set.seed(34)
#' ## Simulates time-series data from the deterministic logistic map
#' ## with a chaotic behaviour.
#' ## ts <- logistic.sim(a=4, s=0, n=1000)
#' ##
#' ## Simulates time-series data from the deterministic logistic map
#' ## with a non-chaotic behaviour.
#' ## ts <- logistic.sim(a=3.2, s=0, n=1000)
#' @author Julio E. Sandubete, Lorenzo Escot
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @importFrom stats as.ts
#' @export logistic.sim
logistic.sim <- function(a = 4, s = 0, x0 = runif(1, 0, 1), n = 1000, n.start = 50) {

  # Settings
  innov <- rnorm(n)
  start.innov <- rnorm(n.start)
  e <- c(start.innov[1L:n.start], innov[1L:n])
  ntot <- length(e)
  x <- double(ntot)
  x[1] <- x0

  # Simulates time-series data from the Logistic map
  if (s != 0) {
    for (i in 2:ntot) {
      x[i] <- a * x[i - 1] * (1 - x[i - 1])
    }
    x <- x + s * e
  } else {
    for (i in 2:ntot) {
      x[i] <- a * x[i - 1] * (1 - x[i - 1])
    }
  }
  if (n.start > 0) {
    x <- x[-(1L:n.start)]
  }

  # Output
  return(ts(x))
}
