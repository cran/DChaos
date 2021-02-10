################################################################################
############ DChaos: An R Package for Chaotic Time Series Analysis #############
################################################################################
# FUNCTION:         DESCRIPTION:
#   netfit            Fits a neural network to an unknown function
################################################################################
#' Fits any standard feedforward neural net model from time-series data
#' @name netfit
#' @aliases netfit
#' @description
#' This function fits any standard feedforward neural net model from time-series data.
#' @param serie a \code{vector}, a time-series object \code{ts} or \code{xts}, a \code{data.frame}, a \code{data.table} or a \code{matrix} depending on the method selected in \code{timelapse}.
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
#' @return This function returns several objects considering the parameter set selected by the user. The best-fitted feed-forward single hidden layer neural net model is saved. It also contains some useful information about the best set of weights found, the fitted values, the residuals obtained or the best embedding parameters set chosen. The best 10 models are displayed on the console. The first column is the neural net number, the second column is the embedding dimension, the third column is the lag or reconstruction delay considered, the fourth column is the number of neurones (or nodes) in the single hidden layer and the fifth column is the Bayesian Information Criterion (BIC) value corresponding to that neural net. Notice that the neural net models are sorted from lowest to highest BIC values.
#' @note The process of adjustment to a neural net model often suffers from being trapped in local optima and different initialization strategies should be taken into account. For this reason the function \code{w0.net} have been implemented. This function estimates previously the initial parameter vector of the neural net model being able to set the maximum number of iterations that the user wants to obtain setting \code{w0maxit}. In addition, by default the neural network estimation is initialized with a fixed seed denoted by \code{seed.t=TRUE} with a value equal to \code{seed=56666459}. The R user can let the seed be fixed either randomly by \code{seed.t=FALSE} or even fix other value of the seed to be able to replicate the results obtained.
#' @references Ripley, B.D. 1996 Pattern Recognition and Neural Networks. Cambridge.
#' @references Venables, W.N., Ripley, B.D. 2002 Modern Applied Statistics with S. Fourth edition. Springer.
#' @references Hornik, K., Stinchcombe, M., White, H. 1989 Multilayer feedforward networks are universal approximators. Neural Networks 2(5):359-366.
#' @examples
#' ## set.seed(34)
#' ## Simulates time-series data from the Logistic map with chaos
#' ## ts        <- DChaos::logistic.sim(n=1000, a=4)
#' ## show(head(ts, 5))
#'
#' ## Provides the best-fitted neural network models for certain parameter set
#' ## model    <- DChaos::netfit(ts, m=1:4, lag=1:3, timelapse="FIXED", h=2:10)
#' ## summary(model)
#' @author Julio E. Sandubete, Lorenzo Escot
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom nnet nnet
#' @importFrom stats residuals
#' @export netfit
netfit <- function(serie, m = 1:4, lag = 1:1, timelapse = c("FIXED", "VARIABLE"), h = 2:10, w0maxit = 100, wtsmaxit = 1e6, pre.white = TRUE, trace = 1, seed.t = TRUE, seed = 56666459, ...) {

  # Checks
  if (is.null(serie)) {
    stop("'serie' should be a vector, a time-series object ts or xts, a data.frame, a data.table or a matrix depending on the method selected in 'timelapse'")
  }
  if (min(m) < 1) {
    stop("wrong value for the embedding dimension")
  }
  if (min(lag) < 1) {
    stop("wrong value for the reconstruction delay")
  }
  if (min(h) < 2) {
    stop("wrong value for the number of neurones (or nodes) in the single hidden layer")
  }
  if (is.null(timelapse)) {
    stop("'timelapse' should be 'FIXED' or 'VARIABLE'")
  }

  # Settings
  BIC <- c()
  PAR <- expand.grid(m = m, lag = lag, h = h)
  np <- nrow(PAR)
  bar <- utils::txtProgressBar(min = 0, max = np, style = 3)

  # Fit neural networks
  for (i in 1:np) {
    xdata <- embedding(serie, m = PAR[i, 1] + 1, lag = PAR[i, 2], timelapse = timelapse)
    xx <- as.data.frame(xdata[, -1])
    yy <- xdata[, 1]
    w0 <- w0.net(xx, yy, m = PAR[i, 1], h = PAR[i, 3], w0maxit = w0maxit, seed.t = seed.t, seed = seed)
    model <- nnet::nnet(xx, yy, size = PAR[i, 3], Wts = w0, linout = TRUE, wtsmaxit = wtsmaxit, trace = trace > 1, ...)
    num <- mean(residuals(model)^2)
    n <- nrow(xx)
    BIC[i] <- log(num) + log(n) / n * (1 + PAR[i, 3] * (PAR[i, 1] + 2))
    utils::setTxtProgressBar(bar, i)
  }
  res <- cbind(PAR, BIC)
  tab <- sort(BIC, index = TRUE)$ix

  # Console output
  if (trace == 1) {
    cat("\n Best models: \n")
    print(res[tab[1:min(10, nrow(res))], ])
  }

  # Fit selected model
  xdata <- embedding(serie, m = PAR[tab[1], 1] + 1, lag = PAR[tab[1], 2], timelapse = timelapse)
  emb.x <- xdata[, -1]
  emb.y <- xdata[, 1]
  w0 <- w0.net(emb.x, emb.y, m = PAR[tab[1], 1], h = PAR[tab[1], 3], w0maxit = w0maxit, seed.t = seed.t, seed = seed)
  model.fit <- nnet::nnet(emb.x, emb.y, size = PAR[tab[1], 3], Wts = w0, linout = TRUE, wtsmaxit = wtsmaxit, trace = trace > 1, ...)
  model.fit <- c(model.fit, emb.x = list(emb.x), emb = PAR[tab[1], c(1:3)])

  # Class definition
  class(model.fit) <- "nnet"

  # Output
  return(model.fit)
}
