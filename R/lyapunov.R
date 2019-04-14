################################################################################
###### DChaos: an R package for detecting chaotic signals in time series #######
################################################################################
# FUNCTION:         DESCRIPTION:
#   lyapunov          Estimates both the largest Lyapunov exponent (NORMA-2) and
#                     the Lyapunov exponent spectrum (QR decomposition) taking
#                     into account the full sample and three different methods
#                     of subsampling by blocks.
################################################################################
#' Estimation of the Lyapunov exponent through several methods
#' @name lyapunov
#' @aliases lyapunov
#' @description
#' This function estimates both the largest Lyapunov exponent through the Norma-2 method and the Lyapunov exponent spectrum through the QR decomposition method taking into account the full sample and three different methods of subsampling by blocks considering the argument set selected by the user.
#' @param x a numeric vector, time serie, data frame or matrix depending on the method selected in \code{timelapse}.
#' @param lag a non-negative integer denoting the reconstruction delay (Default 1).
#' @param timelapse a character denoting if you consider that the observations are sampled at uniform time intervals \code{FIXED} or with a variable time-lapse between each observation \code{VARIABLE} (Default \code{FIXED}).
#' @param M0 a non-negative integer denoting a lower bound for the embedding dimension (Default 3).
#' @param M1 a non-negative integer denoting an upper bound for the embedding dimension (Default 10).
#' @param H0 a non-negative integer denoting a lower bound for the number of neurones in the hidden layer (Default 2).
#' @param H1 a non-negative integer denoting an upper bound for the number of neurones in the hidden layer (Default 10).
#' @param I a non-negative integer denoting a number of neural networks iterations (Default 100).
#' @param lyapmethod a character denoting if you want to estimate the largest Lyapunov exponent \code{LLE}, the Lyapunov exponent spectrum \code{SLE} or both \code{ALL} (Default \code{LLE}).
#' @param blocking a character denoting if you consider the full sample \code{FULL}, the non-overlapping sample \code{NOVER}, the equally spaced sample \code{EQS}, the bootstrap sample \code{BOOT} or all of them \code{ALL} (Default \code{FULL}).
#' @param pre.white a character denoting if you want to use as points to evaluate the partial derivatives the delayed vectors filtered by the neural network \code{TRUE} or not \code{FALSE} (Default \code{TRUE}).
#' @param B a non-negative integer denoting the number of bootstrap iterations (Default 100).
#' @param netplot a logical value denoting if you want to draw a plot \code{TRUE} or not \code{FALSE}. If it is \code{TRUE} shows as many graphs as networks have been considered. Each of them represents the network structure by drawing the weights with positive values in black and the weights with negative values in grey. The thickness of the lines represents a greater or lesser value (Default \code{TRUE}).
#' @param doplot a logical value denoting if you want to draw a plot \code{TRUE} or not \code{FALSE}. If it is \code{TRUE} the evolution of the Lyapunov exponent values are represented for the whole period considering the lyapunov methods and the blocking methods selected by the user (Default \code{TRUE}).
#' @details If \code{FIXED} has been selected \code{x} must be a numeric vector or time serie. Otherwise \code{VARIABLE} has to be specified. In this case \code{x} must be a data frame or matrix with two columns. First, the date with the following format \code{YMD H:M:OS3} considering milliseconds e.g., 20190407 00:00:03.347. If you don't consider milliseconds you must put .000 after the seconds. It should be an object of class \code{Factor}. Second, the univariate time serie as a sequence of numerical values.
#' @return A list containing the largest Lyapunov exponent, the Lyapunov exponent spectrum or both for each neural network structure considered by keeping to \code{Lyapunov.net}. The dataset saved by each blocking method are the estimated Lyapunov exponent value, the standard error, the z-test value and the p-value. If the user chooses the non-overlapping sample, the equally spaced sample or the bootstrap sample the mean and median values of the Lyapunov exponent are showed. Also some details about the embedding dimension, the sample size, the block length and the block number are recorded.
#' @references Ellner, S., Gallant, A., McCaffrey, D., Nychka, D. 1991 Convergence rates and data requirements for jacobian-based estimates of lyapunov exponents from data. Physics Letters A 153(6):357-363.
#' @references McCaffrey, D.F., Ellner, S., Gallant, A.R., Nychka, D.W. 1992 Estimating the lyapunov exponent of a chaotic system with nonparametric regression. Journal of the American Statistical Association 87(419):682-695.
#' @references Nychka, D., Ellner, S., Gallant, A.R., McCaffrey, D. 1992 Finding chaos in noisy systems. Journal of the Royal Statistical Society 54(2):399-426.
#' @references Whang, Y.J., Linton, O. 1999 The asymptotic distribution of nonparametric estimates of the lyapunov exponent for stochastic time series. Journal of Econometrics 91(1):1-42.
#' @references Shintani, M., Linton, O. 2004 Nonparametric neural network estimation of Lyapunov exponents and a direct test for chaos. Journal of Econometrics 120(1):1-33.
#' @examples
#' ## We show below an example considering time series from the
#' ## logistic equation. We have estimated the Lyapunov exponent
#' ## spectrum by each blocking method for an embedding dimension (m=4).
#' data<-logistic.ts(u.min=4,u.max=4,B=100,doplot=FALSE)
#' ts<-data$`Logistic 100`$time.serie
#' lyapu<-lyapunov(ts,lag=1,timelapse="FIXED",M0=4,M1=4,H0=2,H1=7,I=10,
#'        lyapmethod="SLE",blocking="ALL",pre.white=TRUE,B=30,netplot=FALSE,
#'        doplot=FALSE)
#' show(lyapu$Lyapunov.net$Spectrum.full$Exponent)
#' show(lyapu$Lyapunov.net$Spectrum.nonoverlap$Exponent)
#' show(lyapu$Lyapunov.net$Spectrum.equally$Exponent)
#' show(lyapu$Lyapunov.net$Spectrum.bootstrap$Exponent)
#' @author Julio E. Sandubete, Lorenzo Escot
#' @seealso \code{\link{embedding}}, \code{\link{jacobi}}, \code{\link{lyapunov.max}}, \code{\link{lyapunov.spec}}
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
lyapunov<-function(x,lag=1,timelapse="FIXED",M0=3,M1=10,H0=2,H1=10,I=100,lyapmethod=c("LLE","SLE","ALL"),blocking=c("FULL","NOVER","EQS","BOOT","ALL"),pre.white=TRUE,B=100,netplot=TRUE,doplot=TRUE){

  # Checks
  if (is.null(x)){stop("'x' should be a numeric vector, time serie, data frame or matrix depending on the method selected in 'timelapse'")}
  if (lag < 1){stop("wrong value for the reconstruction delay")}
  if (is.null(timelapse)){stop("'timelapse' should be 'FIXED' or 'VARIABLE'")}
  if ((M0 < 3) | (M1 < M0) | M1 > length(x)){stop("wrong value for the embedding dimension")}
  if ((H0 < 2) | (H1 < H0) | H1 > length(x)){stop("wrong value for the number of neurones in the hidden layer")}
  if (M0 == M1){
    net.fit<-matrix(nrow = 1,ncol = 3,dimnames = list(NULL,c("Embedding","Neurones","Bayesian IC")))
  } else {
    net.fit<-matrix(nrow = M1-(M0-1),ncol = 3,dimnames = list(NULL,c("Embedding","Neurones","Bayesian IC")))
  }
  if (I<1){stop("wrong value of neural networks iterations")}
  if (is.null(lyapmethod)){stop("'lyapmethod' should be 'LLE', 'SLE' or 'ALL'")}
  if (is.null(blocking)){stop("'blocking' should be 'FULL', 'NOVER', 'EQS', 'BOOT' or 'ALL'")}
  if (B<1){stop("wrong value of bootstrap iterations")}

  # Settings
  LE.list<-vector("list",nrow(net.fit))
  deriv<-jacobi(x,lag=lag,timelapse=timelapse,M0=M0,M1=M1,H0=H0,H1=H1,I=I,pre.white=pre.white,doplot=netplot)
  for (i in 1:nrow(net.fit)){

    x.lyapunov<-deriv[[1+i]]
    lyapmethod=match.arg(lyapmethod)
    blocking=match.arg(blocking)

    # Estimates the largest Lyapunov Exponent by each blocking method
    if(lyapmethod=="LLE" && blocking=="FULL")
      LE=lyapunov.max(x.lyapunov,blocking=blocking,doplot=doplot)

    if(lyapmethod=="LLE" && blocking=="NOVER")
      LE=lyapunov.max(x.lyapunov,blocking=blocking,doplot=doplot)

    if(lyapmethod=="LLE" && blocking=="EQS")
      LE=lyapunov.max(x.lyapunov,blocking=blocking,doplot=doplot)

    if(lyapmethod=="LLE" && blocking=="BOOT")
      LE=lyapunov.max(x.lyapunov,blocking=blocking,B=B,doplot=doplot)

    if(lyapmethod=="LLE" && blocking=="ALL"){
      LE<-vector("list",4)
      names(LE)<-c("Largest.full","Largest.nonoverlap","Largest.equally","Largest.bootstrap")
      LE[[1]]=lyapunov.max(x.lyapunov,blocking="FULL",doplot=doplot)
      LE[[2]]=lyapunov.max(x.lyapunov,blocking="NOVER",doplot=doplot)
      LE[[3]]=lyapunov.max(x.lyapunov,blocking="EQS",doplot=doplot)
      LE[[4]]=lyapunov.max(x.lyapunov,blocking="BOOT",B=B,doplot=doplot)
    }

    # Estimates the Lyapunov exponent spectrum by each blocking method
    if(lyapmethod=="SLE" && blocking=="FULL")
      LE=lyapunov.spec(x.lyapunov,blocking=blocking,doplot=doplot)

    if(lyapmethod=="SLE" && blocking=="NOVER")
      LE=lyapunov.spec(x.lyapunov,blocking=blocking,doplot=doplot)

    if(lyapmethod=="SLE" && blocking=="EQS")
      LE=lyapunov.spec(x.lyapunov,blocking=blocking,doplot=doplot)

    if(lyapmethod=="SLE" && blocking=="BOOT")
      LE=lyapunov.spec(x.lyapunov,blocking=blocking,B=B,doplot=doplot)

    if(lyapmethod=="SLE" && blocking=="ALL"){
      LE<-vector("list",4)
      names(LE)<-c("Spectrum.full","Spectrum.nonoverlap","Spectrum.equally","Spectrum.bootstrap")
      LE[[1]]=lyapunov.spec(x.lyapunov,blocking="FULL",doplot=doplot)
      LE[[2]]=lyapunov.spec(x.lyapunov,blocking="NOVER",doplot=doplot)
      LE[[3]]=lyapunov.spec(x.lyapunov,blocking="EQS",doplot=doplot)
      LE[[4]]=lyapunov.spec(x.lyapunov,blocking="BOOT",B=B,doplot=doplot)
    }

    # Estimates both the largest Lyapunov exponent and the Lyapunov exponent spectrum by each blocking method
    if(lyapmethod=="ALL" && blocking=="FULL"){
      LE<-vector("list",2)
      names(LE)<-c("Largest.full","Spectrum.full")
      LE[[1]]=lyapunov.max(x.lyapunov,blocking=blocking,doplot=doplot)
      LE[[2]]=lyapunov.spec(x.lyapunov,blocking=blocking,doplot=doplot)
    }

    if(lyapmethod=="ALL" && blocking=="NOVER"){
      LE<-vector("list",2)
      names(LE)<-c("Largest.nonoverlap","Spectrum.nonoverlap")
      LE[[1]]=lyapunov.max(x.lyapunov,blocking=blocking,doplot=doplot)
      LE[[2]]=lyapunov.spec(x.lyapunov,blocking=blocking,doplot=doplot)
    }

    if(lyapmethod=="ALL" && blocking=="EQS"){
      LE<-vector("list",2)
      names(LE)<-c("Largest.equally","Spectrum.equally")
      LE[[1]]=lyapunov.max(x.lyapunov,blocking=blocking,doplot=doplot)
      LE[[2]]=lyapunov.spec(x.lyapunov,blocking=blocking,doplot=doplot)
    }

    if(lyapmethod=="ALL" && blocking=="BOOT"){
      LE<-vector("list",2)
      names(LE)<-c("Largest.bootstrap","Spectrum.bootstrap")
      LE[[1]]=lyapunov.max(x.lyapunov,blocking=blocking,doplot=doplot)
      LE[[2]]=lyapunov.spec(x.lyapunov,blocking=blocking,B=B,doplot=doplot)
    }

    if(lyapmethod=="ALL" && blocking=="ALL"){
      LE<-vector("list",8)
      names(LE)<-c("Largest.full","Largest.nonoverlap","Largest.equally","Largest.bootstrap","Spectrum.full","Spectrum.nonoverlap","Spectrum.equally","Spectrum.bootstrap")
      LE[[1]]=lyapunov.max(x.lyapunov,blocking="FULL",doplot=doplot)
      LE[[2]]=lyapunov.max(x.lyapunov,blocking="NOVER",doplot=doplot)
      LE[[3]]=lyapunov.max(x.lyapunov,blocking="EQS",doplot=doplot)
      LE[[4]]=lyapunov.max(x.lyapunov,blocking="BOOT",B=B,doplot=doplot)
      LE[[5]]=lyapunov.spec(x.lyapunov,blocking="FULL",doplot=doplot)
      LE[[6]]=lyapunov.spec(x.lyapunov,blocking="NOVER",doplot=doplot)
      LE[[7]]=lyapunov.spec(x.lyapunov,blocking="EQS",doplot=doplot)
      LE[[8]]=lyapunov.spec(x.lyapunov,blocking="BOOT",B=B,doplot=doplot)
    }
    LE.list[[i]]<-LE
  }

  # Output
  return(c(list(Network.set=deriv[[1]]),Lyapunov.net=LE.list))
}
