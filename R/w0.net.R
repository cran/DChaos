################################################################################
############ DChaos: An R Package for Chaotic Time Series Analysis #############
################################################################################
# FUNCTION:         DESCRIPTION:
#   w0.net            Estimates the initial parameter vector of the neural net model
################################################################################
#' Estimates the initial parameter vector of the neural net model
#' @name w0.net
#' @aliases w0.net
#' @description
#' This function estimates the initial parameter vector of the neural net model.
#' @param x a \code{matrix} or a \code{data.frame} denoting the explanatory variables.
#' @param y a \code{vector}, a \code{matrix} or a \code{data.frame} denoting the response variable.
#' @param m a non-negative integer denoting the embedding dimension (Default 2).
#' @param h a non-negative integer denoting the number of neurones (or nodes) in the single hidden layer (Default 2).
#' @param rangx a non-negative integer denoting the range of the explanatory variables (Default 1/max(abs(x)).
#' @param w0maxit a non-negative integer denoting the maximum iterations to estimate the initial parameter vector of the neural net models (Default 100).
#' @param seed.t a logical value denoting if the user wants to fix the seed \code{TRUE} or not \code{FALSE} (Default TRUE).
#' @param seed a non-negative integer denoting the value of the seed selected if \code{seed.t = TRUE} (Default 56666459).
#' @return The optimal initial parameter vector of the neural net model considering the argument set selected by the user.
#' @note The process of adjustment to a neural network often suffers from being trapped in local optima and different initialization strategies should be taken into account. For this reason the function \code{w0.net} have been implemented. This function estimates previously the initial parameter vector of the neural net model being able to set the maximum number of iterations that the user wants to obtain setting \code{w0maxit}. In addition, by default the neural network estimation is initialized with a fixed seed denoted by \code{seed.t=TRUE} with a value equal to \code{seed=56666459}. The R user can let the seed be fixed either randomly by \code{seed.t=FALSE} or even fix other value of the seed to be able to replicate the results obtained.
#' @references Ripley, B.D. 1996 Pattern Recognition and Neural Networks. Cambridge.
#' @references Venables, W.N., Ripley, B.D. 2002 Modern Applied Statistics with S. Fourth edition. Springer.
#' @references Hornik, K., Stinchcombe, M., White, H. 1989 Multilayer feedforward networks are universal approximators. Neural Networks 2(5):359-366.
#' @author Julio E. Sandubete, Lorenzo Escot
#' @importFrom nnet nnet
#' @export w0.net
w0.net        <- function(x, y, m=2, h=2, rangx=1/max(abs(x)), w0maxit=100, seed.t=TRUE, seed=56666459){

  # Checks
  if (is.null(x)){stop("'x' should be a numeric vector, time serie, data frame or matrix denoting the explanatory variables")}
  if (is.null(y)){stop("'y' should be a numeric vector, time serie, data frame or matrix denoting the response variable")}
  if (m < 1){stop("wrong value for the embedding dimension")}
  if (h < 2){stop("wrong value for the number of neurones (or nodes) in the single hidden layer")}

  # Settings
  Wts0        <- NULL
  ECMWts0     <- NULL
  if(seed.t==T){
    set.seed(seed)
  } else {
    set.seed()
  }

  # Estimates the initial parameter vector of the neural net model
  for (i in 1:w0maxit){
    W0        <- runif(h*(m+1)+h+1, min=-rangx, max=rangx)
    net.nn    <- nnet::nnet(x,y,size=h,Wts=W0,maxit=0,linout = T,trace=F)
    if(i==1){
      ECMWts0 <- net.nn$value
      Wts0    <- W0
    }
    if(net.nn$value<ECMWts0){
      ECMWts0 <- net.nn$value
      Wts0    <- W0
    }
  }

  # Output
  return(Wts0)
}
