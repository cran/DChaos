################################################################################
############ DChaos: An R Package for Chaotic Time Series Analysis #############
################################################################################
# FUNCTION:             DESCRIPTION:
#   summary.lyapunov      Summarizing estimated Lyapunov exponents
################################################################################
#' Summary method for a lyapunov object
#' @name summary.lyapunov
#' @aliases summary.lyapunov
#' @description
#' summary method for class "lyapunov".
#' @usage \method{summary}{lyapunov}(object, ...)
#' @param object an object of class \code{"lyapunov"} provided by \code{\link{lyapunov.max}}, \code{\link{lyapunov.spec}} or \code{\link{lyapunov}} functions.
#' @param ... further arguments passed to or from other methods.
#' @return This function \code{summary.lyapunov} computes and returns a list of summary statistics of the results given in a \code{lyapunov} object using the components (list elements) from its argument.
#' @examples
#' ## set.seed(34)
#' ## Simulates time-series data from the Logistic map with chaos
#' ## ts        <- DChaos::logistic.sim(n=1000, a=4)
#' ## show(head(ts, 5))
#'
#' ## Summary method for a lyapunov object (only 1 method)
#' ## jacobian <- DChaos::jacobian.net(data=ts, m=3:3, lag=1:1, timelapse="FIXED", h=2:10)
#' ## exponent <- DChaos::lyapunov.spec(data=jacobian, blocking="BOOT", B=100, doplot=FALSE)
#' ## summary(exponent)
#'
#' ## Summary method for a lyapunov object (> 1 method)
#' ## exponent <- DChaos::lyapunov(ts, m=3:3, lag=1:1, timelapse="FIXED", h=2:10, w0maxit=100,
#' ##                     wtsmaxit=1e6, pre.white=TRUE, lyapmethod="SLE", blocking="ALL",
#' ##                     B=100, trace=1, seed.t=TRUE, seed=56666459, doplot=FALSE))
#' ## summmary(exponent)
#' @author Julio E. Sandubete, Lorenzo Escot
#' @export
summary.lyapunov <- function(object, ...) {
  if (object$nprint == 0) {
    cat("Call:\n")
    cat(object$estimator, "\n")
    cat("\nCoefficients:\n")
    if (object$procedure == "QR decomposition by full sample method" | object$procedure == "Norma-2 by full sample method") {
      print(object$exponent)
    } else {
      print(object$exponent.median)
    }
    cat("---\n")
    cat("Procedure:", object$procedure, "\n")
    cat("Embedding dimension: ", object$emb.m, ", ", "Time-delay: ", object$emb.lag, ", ", "No. hidden units: ", object$emb.h, sep = "")
    cat("\nSample size: ", object$sample, ", ", "Block length: ", object$block.length, ", ", "No. blocks: ", object$no.block, sep = "")
  } else {
    cat("Call:\n")
    cat(unlist(object[[21]][1]), "\n")
    cat("\nCoefficients:\n")
    if (unlist(object[[21]][[2]][1]) == "QR decomposition by full sample method" | unlist(object[[21]][[2]][1]) == "Norma-2 by full sample method") {
      print(object[[21]][[3]])
      cat("---\n")
      cat("Procedure:", unlist(object[[21]][2]), "\n")
      cat("Embedding dimension: ", object$emb.m, ", ", "Time-delay: ", object$emb.lag, ", ", "No. hidden units: ", object$emb.h, sep = "")
      cat("\nSample size: ", unlist(object[[21]][4]), ", ", "Block length: ", unlist(object[[21]][5]), ", ", "No. blocks: ", unlist(object[[21]][6]), sep = "")
      cat("... only the first method is shown (see lyapunov object)\n")
    } else {
      print(object[[21]][[4]])
      cat("---\n")
      cat("Procedure:", unlist(object[[21]][2]), "\n")
      cat("Embedding dimension: ", object$emb.m, ", ", "Time-delay: ", object$emb.lag, ", ", "No. hidden units: ", object$emb.h, sep = "")
      cat("\nSample size: ", unlist(object[[21]][5]), ", ", "Block length: ", unlist(object[[21]][6]), ", ", "No. blocks: ", unlist(object[[21]][7]), sep = "")
      cat("... only the first method is shown (see lyapunov object)\n")
    }
  }
}
