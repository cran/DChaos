################################################################################
############ DChaos: An R Package for Chaotic Time Series Analysis #############
################################################################################
# FUNCTION:         DESCRIPTION:
#   lyapunov          Estimates the Lyapunov exponent through several methods
################################################################################
#' Estimates the Lyapunov exponent through several methods
#' @name lyapunov
#' @aliases lyapunov
#' @description
#' This is an all-in-one function. It provides, at the same time, the delayed-coordinate embedding vector (\code{embedding}), estimates the best neural net model (\code{netfit}), calculates the partial derivatives directly from the chosen neural network model (\code{jacobian.net}). Finally, this function estimates both the largest Lyapunov exponent through the Norma-2 procedure (\code{lyapunov.max}) and the Lyapunov exponent spectrum through the QR decomposition procedure (\code{lyapunov.spec}) taking into account the full sample and three different methods of subsampling by blocks.
#' @param data a \code{vector}, a time-series object \code{ts} or \code{xts}, a \code{data.frame}, a \code{data.table} or a \code{matrix} depending on the method selected in \code{timelapse}.
#' @param m a non-negative integer denoting a lower and upper bound for the embedding dimension (Default 1:4).
#' @param lag a non-negative integer denoting a lower and upper bound for the the reconstruction delay (Default 1:1).
#' @param timelapse a character denoting if the time-series data are sampled at uniform time-frequency e.g., 1-month, 1-day, 1-hour, 30-min, 5-min, 1-min and so on \code{FIXED} or non-uniform time-frequency which are not equally spaced in time \code{VARIABLE} (Default \code{FIXED}).
#' @param h a non-negative integer denoting a lower and upper bound for the number of neurones (or nodes) in the single hidden layer (Default 2:10).
#' @param w0maxit a non-negative integer denoting the maximum iterations to estimate the initial parameter vector of the neural net models (Default 100).
#' @param wtsmaxit a non-negative integer denoting the maximum iterations to estimate the weights parameter vector of the neural net models (Default 1e6).
#' @param pre.white a logical value denoting if the user wants to use as points to evaluate the partial derivatives the delayed vectors filtered by the neural net model chosen \code{TRUE} or not \code{FALSE} (Default \code{TRUE}).
#' @param lyapmethod a character denoting the procedure chosen to estimate the Lyapunov exponent. If \code{LLE} has been selected the function will estimate only the largest Lyapunov exponent through the Norma-2 method. If \code{SLE} has been selected the function will estimate the Lyapunov exponent spectrum through the QR decomposition. Otherwise \code{ALL} has to be specified. In this case the function will estimate the Lyapunov exponent considering both procedures (Default \code{SLE}).
#' @param blocking a character denoting the blocking method chosen for figuring out the Lyapunov exponent. Available options are \code{FULL} if the user considers the full sample, \code{NOVER} if the user considers the non-overlapping sample, \code{EQS} if the user considers the equally spaced sample, \code{BOOT} if the user considers the bootstrap sample or \code{ALL} if the user considers all of them (Default \code{BOOT}).
#' @param B a non-negative integer denoting the number of bootstrap iterations (Default 1000).
#' @param trace a binary value denoting if the user wants to print the output on the console \code{1} or not \code{0} (Default 1).
#' @param seed.t a logical value denoting if the user wants to fix the seed \code{TRUE} or not \code{FALSE} (Default TRUE).
#' @param seed a non-negative integer denoting the value of the seed selected if \code{seed.t = TRUE} (Default 56666459).
#' @param doplot a logical value denoting if the user wants to draw a plot \code{TRUE} or not \code{FALSE}. If it is \code{TRUE} the evolution of the Lyapunov exponent values are represented for the whole period considering the blocking method chosen by the user. It shows as many graphs as embedding dimensions have been considered (Default \code{TRUE}).
#' @param ... further arguments passed to or from \code{nnet} function.
#' @return This function returns several objects considering the parameter set selected by the user. The largest Lyapunov exponent (Norma-2 procedure) and the Lyapunov exponent spectrum (QR decomposition procedure) by each blocking method are estimated. It also contains some useful information about the estimated jacobian, the best-fitted feed-forward single hidden layer neural net model, the best set of weights found, the fitted values, the residuals obtained, the best embedding parameters set chosen, the sample size or the block length considered by each blocking method. This function provides the standard error, the z test value and the p-value for testing the null hypothesis \eqn{H0: \lambda_k > 0 for k = 1,2,3, \ldots, m}. Reject the null hypothesis ${H_0}$ means lack of chaotic behaviour. That is, the data-generating process does not have a chaotic attractor because of it does not show the property of sensitivity to initial conditions.
#' @note We have considered it appropriate to incorporate a function that unifies the whole process to make it easier and more intuitive for the R users. The DChaos package provides several ways to figure out robustly the neural net estimator of the k-th Lyapunov exponent. Particularly, there are 8 functions (one for each procedure and blocking method) which estimate the Lyapunov exponents consistently. Hence the DChaos package allows the R users to choose between two different procedures to obtain the neural net estimator of the k-th Lyapunov exponent and four ways of subsampling by blocks: full sample, non-overlapping sample, equally spaced sample and bootstrap sample. The blocking methods what they do is to split the time-series data into several blocks by estimating a Lyapunov exponent for each subsample. If the R users choose the non-overlapping sample (\code{blocking = "NOVER"}), the equally spaced sample (\code{blocking = "EQS"}) or the bootstrap sample (\code{blocking = "BOOT"}) the mean and median values of the Lyapunov exponent for each block or subsample are saved. By default we recommend using the median value as it is more robust to the presence of outliers. Notice that the parameter \code{B} will only be considered if the R users choose the bootstrap blocking method.
#' @references Ellner, S., Gallant, A., McCaffrey, D., Nychka, D. 1991 Convergence rates and data requirements for jacobian-based estimates of lyapunov exponents from data. Physics Letters A 153(6):357-363.
#' @references McCaffrey, D.F., Ellner, S., Gallant, A.R., Nychka, D.W. 1992 Estimating the lyapunov exponent of a chaotic system with nonparametric regression. Journal of the American Statistical Association 87(419):682-695.
#' @references Nychka, D., Ellner, S., Gallant, A.R., McCaffrey, D. 1992 Finding chaos in noisy systems. Journal of the Royal Statistical Society 54(2):399-426.
#' @references Whang, Y.J., Linton, O. 1999 The asymptotic distribution of nonparametric estimates of the lyapunov exponent for stochastic time series. Journal of Econometrics 91(1):1-42.
#' @references Shintani, M., Linton, O. 2004 Nonparametric neural network estimation of Lyapunov exponents and a direct test for chaos. Journal of Econometrics 120(1):1-33.
#' @examples
#' ## set.seed(34)
#' ## Simulates time-series data from the Logistic map with chaos
#' ## ts        <- DChaos::logistic.sim(n=1000, a=4)
#' ## show(head(ts, 5))
#'
#' ## Provides the Lyapunov exponent spectrum by the QR decomposition procedure considering the
#' ## bootstrap blocking method directly from the Logistic map with chaos simulated.
#' ## exponent <- DChaos::lyapunov(ts, m=3:3, lag=1:1, timelapse="FIXED", h=2:10, w0maxit=100,
#' ##                     wtsmaxit=1e6, pre.white=TRUE, lyapmethod="SLE", blocking="ALL",
#' ##                     B=100, trace=1, seed.t=TRUE, seed=56666459, doplot=FALSE))
#' ## summmary(exponent)
#' @author Julio E. Sandubete, Lorenzo Escot
#' @seealso \code{\link{lyapunov.max}}, \code{\link{lyapunov.spec}}
#' @importFrom pracma normest
#' @importFrom pracma zeros
#' @importFrom sandwich lrvar
#' @importFrom stats pnorm
#' @importFrom graphics abline
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom stats median
#' @importFrom stats is.ts
#' @export lyapunov
lyapunov <- function(data, m = 1:4, lag = 1:1, timelapse = c("FIXED", "VARIABLE"), h = 2:10, w0maxit = 100, wtsmaxit = 1e6, pre.white = TRUE, lyapmethod = c("SLE", "LLE", "ALL"), blocking = c("BOOT", "NOVER", "EQS", "FULL", "ALL"), B = 1000, trace = 1, seed.t = TRUE, seed = 56666459, doplot = TRUE, ...) {

  # Checks
  if (is.null(data)) {
    stop("'data' should be a vector, a time-series object ts or xts, a data.frame, a data.table or a matrix depending on the method selected in 'timelapse'")
  }
  if (min(m) < 1) {
    stop("wrong value for the embedding dimension")
  }
  if (min(lag) < 1) {
    stop("wrong value for the reconstruction delay")
  }
  if (min(h) < 2) {
    stop("wrong value for the number of neurones in the hidden layer")
  }
  if (is.null(timelapse)) {
    stop("'timelapse' should be 'FIXED' or 'VARIABLE'")
  }
  if (w0maxit < 1) {
    stop("wrong value of neural networks iterations")
  }
  if (B < 1) {
    stop("wrong value of bootstrap iterations")
  }
  if (is.null(lyapmethod)) {
    stop("'lyapmethod' should be 'SLE', 'LLE', or 'ALL'")
  }
  if (is.null(blocking)) {
    stop("'blocking' should be 'BOOT', 'NOVER', 'EQS', 'FULL' or 'ALL'")
  }
  if (B < 1) {
    stop("wrong value of bootstrap iterations")
  }

  # Settings
  jacobian <- jacobian.net(data = data, m = m, lag = lag, timelapse = timelapse, h = h, w0maxit = w0maxit, wtsmaxit = wtsmaxit, pre.white = pre.white, trace = trace, seed.t = seed.t, seed = seed, ...)
  lyapmethod <- match.arg(lyapmethod)
  blocking <- match.arg(blocking)

  # Estimates the Lyapunov exponent spectrum by each blocking method
  if (lyapmethod == "SLE" && blocking == "BOOT") {
    LE <- lyapunov.spec(data = jacobian, blocking = blocking, B = B, doplot = doplot)
  }

  if (lyapmethod == "SLE" && blocking == "NOVER") {
    LE <- lyapunov.spec(data = jacobian, blocking = blocking, doplot = doplot)
  }

  if (lyapmethod == "SLE" && blocking == "EQS") {
    LE <- lyapunov.spec(data = jacobian, blocking = blocking, doplot = doplot)
  }

  if (lyapmethod == "SLE" && blocking == "FULL") {
    LE <- lyapunov.spec(data = jacobian, blocking = blocking, doplot = doplot)
  }

  if (lyapmethod == "SLE" && blocking == "ALL") {
    exponent.boot <- lyapunov.spec(data = jacobian, blocking = "BOOT", B = B, doplot = doplot)
    exponent.nover <- lyapunov.spec(data = jacobian, blocking = "NOVER", doplot = doplot)
    exponent.eqs <- lyapunov.spec(data = jacobian, blocking = "EQS", doplot = doplot)
    exponent.full <- lyapunov.spec(data = jacobian, blocking = "FULL", doplot = doplot)
    LE <- c(jacobian,
      exponent.boot = list(exponent.boot[c(21:27)]), exponent.nover = list(exponent.nover[c(21:27)]),
      exponent.eqs = list(exponent.eqs[c(21:27)]), exponent.full = list(exponent.full[c(21:26)]), nprint = 1
    )
  }

  # Estimates the largest Lyapunov Exponent by each blocking method
  if (lyapmethod == "LLE" && blocking == "BOOT") {
    LE <- lyapunov.max(data = jacobian, blocking = blocking, B = B, doplot = doplot)
  }

  if (lyapmethod == "LLE" && blocking == "NOVER") {
    LE <- lyapunov.max(data = jacobian, blocking = blocking, doplot = doplot)
  }

  if (lyapmethod == "LLE" && blocking == "EQS") {
    LE <- lyapunov.max(data = jacobian, blocking = blocking, doplot = doplot)
  }

  if (lyapmethod == "LLE" && blocking == "FULL") {
    LE <- lyapunov.max(data = jacobian, blocking = blocking, doplot = doplot)
  }

  if (lyapmethod == "LLE" && blocking == "ALL") {
    exponent.boot <- lyapunov.max(data = jacobian, blocking = "BOOT", B = B, doplot = doplot)
    exponent.nover <- lyapunov.max(data = jacobian, blocking = "NOVER", doplot = doplot)
    exponent.eqs <- lyapunov.max(data = jacobian, blocking = "EQS", doplot = doplot)
    exponent.full <- lyapunov.max(data = jacobian, blocking = "FULL", doplot = doplot)
    LE <- c(jacobian,
      exponent.boot = list(exponent.boot[c(21:27)]), exponent.nover = list(exponent.nover[c(21:27)]),
      exponent.eqs = list(exponent.eqs[c(21:27)]), exponent.full = list(exponent.full[c(21:26)]), nprint = 1
    )
  }

  # Estimates both estimators by each blocking method
  if (lyapmethod == "ALL" && blocking == "BOOT") {
    exponent.spec <- lyapunov.spec(data = jacobian, blocking = blocking, B = B, doplot = doplot)
    exponent.max <- lyapunov.max(data = jacobian, blocking = blocking, B = B, doplot = doplot)
    LE <- c(jacobian, exponent.spec = list(exponent.spec[c(21:27)]), exponent.max = list(exponent.max[c(21:27)]), nprint = 1)
  }

  if (lyapmethod == "ALL" && blocking == "NOVER") {
    exponent.spec <- lyapunov.spec(data = jacobian, blocking = blocking, B = B, doplot = doplot)
    exponent.max <- lyapunov.max(data = jacobian, blocking = blocking, B = B, doplot = doplot)
    LE <- c(jacobian, exponent.spec = list(exponent.spec[c(21:27)]), exponent.max = list(exponent.max[c(21:27)]), nprint = 1)
  }

  if (lyapmethod == "ALL" && blocking == "EQS") {
    exponent.spec <- lyapunov.spec(data = jacobian, blocking = blocking, B = B, doplot = doplot)
    exponent.max <- lyapunov.max(data = jacobian, blocking = blocking, B = B, doplot = doplot)
    LE <- c(jacobian, exponent.spec = list(exponent.spec[c(21:27)]), exponent.max = list(exponent.max[c(21:27)]), nprint = 1)
  }

  if (lyapmethod == "ALL" && blocking == "FULL") {
    exponent.spec <- lyapunov.spec(data = jacobian, blocking = blocking, B = B, doplot = doplot)
    exponent.max <- lyapunov.max(data = jacobian, blocking = blocking, B = B, doplot = doplot)
    LE <- c(jacobian, exponent.spec = list(exponent.spec[c(21:26)]), exponent.max = list(exponent.max[c(21:26)]), nprint = 1)
  }

  # Estimates both estimators by al blocking methods
  if (lyapmethod == "ALL" && blocking == "ALL") {
    lyapspec.boot <- lyapunov.spec(data = jacobian, blocking = "BOOT", B = B, doplot = doplot)
    lyapspec.nover <- lyapunov.spec(data = jacobian, blocking = "NOVER", doplot = doplot)
    lyapspec.eqs <- lyapunov.spec(data = jacobian, blocking = "EQS", doplot = doplot)
    lyapspec.full <- lyapunov.spec(data = jacobian, blocking = "FULL", doplot = doplot)
    lyapmax.boot <- lyapunov.max(data = jacobian, blocking = "BOOT", B = B, doplot = doplot)
    lyapmax.nover <- lyapunov.max(data = jacobian, blocking = "NOVER", doplot = doplot)
    lyapmax.eqs <- lyapunov.max(data = jacobian, blocking = "EQS", doplot = doplot)
    lyapmax.full <- lyapunov.max(data = jacobian, blocking = "FULL", doplot = doplot)
    LE <- c(jacobian,
      lyapspec.boot = list(lyapspec.boot[c(21:27)]), lyapmax.boot = list(lyapmax.boot[c(21:27)]),
      lyapspec.nover = list(lyapspec.nover[c(21:27)]), lyapmax.nover = list(lyapmax.nover[c(21:27)]),
      lyapspec.eqs = list(lyapspec.boot[c(21:27)]), lyapmax.eqs = list(lyapmax.boot[c(21:27)]),
      lyapspec.full = list(lyapspec.boot[c(21:26)]), lyapmax.full = list(lyapmax.boot[c(21:26)]), nprint = 1
    )
  }

  # Class definition
  class(LE) <- "lyapunov"

  # Output
  return(LE)
}
