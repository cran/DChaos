################################################################################
############ DChaos: An R Package for Chaotic Time Series Analysis #############
################################################################################
# FUNCTION:         DESCRIPTION:
#   jacobian.net      Partial derivatives are calculated from the best-fitted neural net model
################################################################################
#' Computes the partial derivatives from the best-fitted neural net model
#' @name jacobian.net
#' @aliases jacobian.net
#' @description
#' This function computes analytically the partial derivatives from the best-fitted neural net model.
#' @param model a neural network model fitted using the \code{netfit} function.
#' @param data a \code{vector}, a time-series object \code{ts} or \code{xts}, a \code{data.frame}, a \code{data.table} or a \code{matrix} depending on the method selected in \code{timelapse}.
#' @param m a non-negative integer denoting a lower and upper bound for the embedding dimension (Default 1:4).
#' @param lag a non-negative integer denoting a lower and upper bound for the the reconstruction delay (Default 1:1).
#' @param timelapse a character denoting if the time-series data are sampled at uniform time-frequency e.g., 1-month, 1-day, 1-hour, 30-min, 5-min, 1-min and so on \code{FIXED} or non-uniform time-frequency which are not equally spaced in time \code{VARIABLE} (Default \code{FIXED}).
#' @param h a non-negative integer denoting a lower and upper bound for the number of neurones (or nodes) in the single hidden layer (Default 2:10).
#' @param w0maxit a non-negative integer denoting the maximum iterations to estimate the initial parameter vector of the neural net models (Default 100).
#' @param wtsmaxit a non-negative integer denoting the maximum iterations to estimate the weights parameter vector of the neural net models (Default 1e6).
#' @param pre.white a logical value denoting if the user wants to use as points to evaluate the partial derivatives the delayed vectors filtered by the neural net model chosen \code{TRUE} or not \code{FALSE} (Default \code{TRUE}).
#' @param trace a binary value denoting if the user wants to print the output on the console \code{1} or not \code{0} (Default 1).
#' @param seed.t a logical value denoting if the user wants to fix the seed \code{TRUE} or not \code{FALSE} (Default TRUE).
#' @param seed a non-negative integer denoting the value of the seed selected if \code{seed.t = TRUE} (Default 56666459).
#' @param ... further arguments passed to or from \code{nnet} function.
#' @return This function returns several objects considering the parameter set selected by the user. Partial derivatives are calculated analytically from the best-fitted neural net model. It also contains some useful information about the best-fitted feed-forward single hidden layer neural net model saved, the best set of weights found, the fitted values, the residuals obtained or the best embedding parameters set chosen. This function allows the R user uses the data previously obtained from the best-fitted neural network estimated by the \code{netfit} function if \code{model} is not empty. Otherwise \code{data} has to be specified.
#' @note The main reason for using neural network models is not to look for the best predictive model but to estimate a model that captures the non-linear time dependence well enough and, additionally, allows us to obtain in an analytical way (instead of numerical) the jacobian functional of the unknown underlying generator system. The estimation of this jacobian or partial derivatives will later allow us to contrast our hypothesis of chaos estimating the Lyapunov exponents.
#' @references Eckmann, J.P., Ruelle, D. 1985 Ergodic theory of chaos and strange attractors. Rev Mod Phys 57:617–656.
#' @references Gencay, R., Dechert, W.D. 1992 An algorithm for the n lyapunov exponents of an n-dimensional unknown dynamical system. Physica D 59(1):142–157.
#' @references Shintani, M., Linton, O. 2004 Nonparametric neural network estimation of Lyapunov exponents and a direct test for chaos. Journal of Econometrics 120(1):1-33.
#' @examples
#' ## set.seed(34)
#' ## Simulates time-series data from the Logistic map with chaos
#' ## ts        <- DChaos::logistic.sim(n=1000, a=4)
#' ## show(head(ts, 5))
#'
#' ## Computes analytically the partial derivatives from the best-fitted neural net model
#' ## showed in the netfit example
#' ## model    <- DChaos::netfit(ts, m=1:4, lag=1:3, timelapse="FIXED", h=2:10)
#' ## jacobian <- DChaos::jacobian.net(model=model)
#' ## summary(jacobian)
#'
#' ## Partial derivatives are calculated analytically without setting previously any neural net model
#' ## jacobian <- DChaos::jacobian.net(data=ts, m=3:3, lag=1:1, timelapse="FIXED", h=2:10)
#' ## summary(jacobian)
#' @author Julio E. Sandubete, Lorenzo Escot
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom nnet nnet
#' @export jacobian.net
jacobian.net <- function(model, data, m = 1:4, lag = 1:1, timelapse = c("FIXED", "VARIABLE"), h = 2:10, w0maxit = 100, wtsmaxit = 1e6, pre.white = TRUE, trace = 1, seed.t = TRUE, seed = 56666459, ...) {

  # Checks
  if (missing(model)) {
    model <- netfit(serie = data, m = m, lag = lag, timelapse = timelapse, h = h, w0maxit = w0maxit, wtsmaxit = wtsmaxit, trace = trace, seed.t = seed.t, seed = seed, ...)
  }

  # Settings
  m <- model$n[1]
  h <- model$n[2]
  w <- model$wts
  xo.net <- as.matrix(model$emb.x)
  n <- nrow(xo.net)
  z <- matrix(nrow = n, ncol = h)
  dphi <- matrix(nrow = n, ncol = h)
  zout <- matrix(nrow = n, ncol = h)
  dzout <- matrix(nrow = n, ncol = m)
  wout <- c()
  wz <- c()

  # Calculation and evaluation of partial derivatives (analytically)
  for (i in 1:h) {
    z[, i] <- w[i + (i - 1) * m] + xo.net %*% w[(i + (i - 1) * m + 1):((i + (i - 1) * m) + m)]
    dphi[, i] <- exp(z[, i]) / (1 + exp(z[, i]))^2
    wout[i] <- w[h + 1 + i + h * m]
    zout[, i] <- dphi[, i] * wout[i]
  }
  for (j in 1:m) {
    for (s in 1:h) {
      wz[s] <- w[(s - 1) * m + (s - 1) + (j + 1)]
      wzout <- cbind(wz)
    }
    dzout[, j] <- zout %*% wzout
  }
  dzout <- as.data.frame(dzout[complete.cases(dzout), ])
  colnames(dzout) <- paste0("dx", 1:m)
  model <- c(model, jacobian = list(dzout))

  # Class definition
  class(model) <- "nnet"

  # Output
  return(model)
}
