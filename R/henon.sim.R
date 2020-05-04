################################################################################
############ DChaos: An R Package for Chaotic Time Series Analysis #############
################################################################################
# FUNCTION:         DESCRIPTION:
#   henon.sim         Simulates time series from the Henon map
################################################################################
#' Simulates time-series data from the Henon map
#' @name henon.sim
#' @aliases henon.sim
#' @description
#' This function simulates time-series data from the Henon map considering the parameter set selected by the user. The initial condition is a random number between -0.5 and 0.5. Some initial conditions may lead to an unstable system that will tend to infinity.
#' @param a a non-negative integer denoting the value of parameter \code{a} (Default 1.4).
#' @param b a non-negative integer denoting the value of parameter \code{b} (Default 0.3).
#' @param s a non-negative integer denoting the variance value of the error term. If \eqn{s=0} gives the standard deterministic map (Default 0).
#' @param x0 a non-negative integer denoting the initial condition of x-coordinate (Default random number between -0.5 and 0.5).
#' @param y0 a non-negative integer denoting the initial condition of y-coordinate (Default random number between -0.5 and 0.5).
#' @param n a non-negative integer denoting the length (Default 1000).
#' @param n.start a non-negative integer denoting the number of observations that will be discarded to ensure that the values are in the attractor (Default 50).
#' @return A time-series data object generated from the Henon map with or without an additive measurement noise term. This dataset could be useful for researchers interested in the field of chaotic dynamic systems and non-linear time series analysis and professors (and students) who teach (learn) courses related to those topics.
#' @note This function provides also noisy time-series data from the deterministic henon map adding an additive measurement noise term if \eqn{s>0}. We have added to each time-series data a normal multinomial error term denoted by \eqn{{\varepsilon _t} \sim N\left( {0,s} \right)} with different variance values (\eqn{s}). In this sense we have considered it appropriate to add a measurement noise term because most real-world observed time-series data are usually noise-contaminated signals, characterised by an erratic and persistent volatility in certain periods and there is almost always a source of noise linked to measurement errors in real-world datasets.
#' @references Hénon, M. 1976 A two-dimensional mapping with a strange attractor. Communications in Mathematical Physics 50(1):69-77.
#' @examples
#' ## set.seed(34)
#' ## Simulates time-series data from the deterministic henon map
#' ## with a chaotic behaviour.
#' ts <- henon.sim(a=1.4, b=0.3, s=0, n=1000)
#' ##
#' ## Simulates time-series data from the deterministic henon map
#' ## with a non-chaotic behaviour.
#' ts <- henon.sim(a=1.2, b=0.1, s=0, n=1000)
#' @author Julio E. Sandubete, Lorenzo Escot
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @importFrom stats as.ts
#' @export henon.sim
henon.sim     <- function(a=1.4, b=0.3, s=0, x0=runif(1,-0.5,0.5), y0=runif(1,-0.5,0.5), n=1000, n.start = 50){

  # Settings
  innov         <- rnorm(n)
  start.innov   <- rnorm(n.start)
  e             <- c(start.innov[1L:n.start], innov[1L:n])
  ntot          <- length(e)
  x             <- double(ntot)
  y             <- double(ntot)
  x[1]          <- x0
  y[1]          <- y0

  # Simulates time-series data from the henon map
  if(s != 0) {
    for(i in 2:ntot) {
      x[i] <- 1-a*x[i-1]^2+b*y[i-1]
      y[i] <- x[i-1]
    }
    x <- x + s*e
    y <- y + s*e
  } else {
    for(i in 2:ntot) {
      x[i] <- 1-a*x[i-1]^2+b*y[i-1]
      y[i] <- x[i-1]
    }
  }
  if (n.start > 0){
    x <- x[-(1L:n.start)]
    y <- y[-(1L:n.start)]
  }

  # Output
  return(ts(cbind(x,y)))
}


