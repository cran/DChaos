################################################################################
###### DChaos: an R package for detecting chaotic signals in time series #######
################################################################################
# FUNCTION:         DESCRIPTION:
#   jacobi            Applicates the Jacobian method by a fit through neural
#                     networks
################################################################################
#' Application of Jacobian method by a fit through neural networks
#' @name jacobi
#' @aliases jacobi
#' @description
#' This function estimates the partial derivatives of the jacobian by a fit through a single-hidden-layer neural network considering the argument set selected by the user.
#' @param x a numeric vector, time serie, data frame or matrix depending on the method selected in \code{timelapse}.
#' @param lag a non-negative integer denoting the reconstruction delay (Default 1).
#' @param timelapse a character denoting if you consider that the observations are sampled at uniform time intervals \code{FIXED} or with a variable time-lapse between each observation \code{VARIABLE} (Default \code{FIXED}).
#' @param M0 a non-negative integer denoting a lower bound for the embedding dimension (Default 3).
#' @param M1 a non-negative integer denoting an upper bound for the embedding dimension (Default 10).
#' @param H0 a non-negative integer denoting a lower bound for the number of neurones in the hidden layer (Default 2).
#' @param H1 a non-negative integer denoting an upper bound for the number of neurones in the hidden layer (Default 10).
#' @param I a non-negative integer denoting a number of neural networks iterations (Default 100).
#' @param pre.white a character denoting if you want to use as points to evaluate the partial derivatives the delayed vectors filtered by the neural network \code{TRUE} or not \code{FALSE} (Default \code{TRUE}).
#' @param doplot a logical value denoting if you want to draw a plot \code{TRUE} or not \code{FALSE}. If it is \code{TRUE} shows as many graphs as networks have been considered. Each of them represents the network structure by drawing the weights with positive values in black and the weights with negative values in grey. The thickness of the lines represents a greater or lesser value (Default \code{TRUE}).
#' @details If \code{FIXED} has been selected \code{x} must be a numeric vector or time serie. Otherwise \code{VARIABLE} has to be specified. In this case \code{x} must be a data frame or matrix with two columns. First, the date with the following format \code{YMD H:M:OS3} considering milliseconds e.g., 20190407 00:00:03.347. If you don't consider milliseconds you must put .000 after the seconds. It should be an object of class \code{Factor}. Second, the univariate time serie as a sequence of numerical values.
#' @return A list with several objects. The first output is a matrix called \code{Network.set}. It contains the networks that have the best fit for each embedding dimension \code{m}. That is, the neural networks that have the minimum bayesian information criterion (BIC) between all possible number of neurones in the hidden layer. Then, the partial derivatives of the jacobian are saved on a data frame for each neural network structure considered by keeping to \code{Jacobian.net}.
#' @references Eckmann, J.P., Ruelle, D. 1985 Ergodic theory of chaos and strange attractors. Reviews of Modern Physics 57:617-656.
#' @references Eckmann, J.P., Kamphorst, S.O., Ruelle, D., Ciliberto, S. 1986 Liapunov exponents from time series. Physical Review A 34:971-979.
#' @references Hornik, K., Stinchcombe, M., White, H. 1989 Multilayer feedforward networks are universal approximators. Neural Networks 2(5):359-366.
#' @references Gencay, R., Dechert, W. 1992 An algorithm for the n lyapunov exponents of an n-dimensional unknown dynamical system. Physica D 59(1):142-157.
#' @references McCaffrey, D.F., Ellner, S., Gallant, A.R., Nychka, D.W. 1992 Estimating the lyapunov exponent of a chaotic system with nonparametric regression. Journal of the American Statistical Association 87(419):682-695.
#' @references Nychka, D., Ellner, S., Gallant, A.R., McCaffrey, D. 1992 Finding chaos in noisy systems. Journal of the Royal Statistical Society 54(2):399-426.
#' @references Kuan, C., Liu, T., Gencay, R. 2004 Netfile 4.01: Feedforward neural networks and Lyapunov exponents estimation. Ball State University.
#' @examples
#' ## We show below an example considering time series from the
#' ## logistic equation. The first objetc is a matrix called
#' ## Network.set. It contains the networks that have the best
#' ## fit for each embedding dimension (3<m<4).
#' data<-logistic.ts(u.min=4,u.max=4,B=100,doplot=FALSE)
#' ts<-data$`Logistic 100`$time.serie
#' jacob<-jacobi(ts,lag=1,timelapse="FIXED",M0=3,M1=4,
#'        H0=3,H1=7,I=10,pre.white=TRUE,doplot=FALSE)
#' show(jacob$Network.set)
#' ## The partial derivatives of the jacobian are saved on a
#' ## data frame for each neural network structure considered
#' ## by keeping to Jacobian.net. The first ten jacobian values
#' ## corresponding to the neural network for m=4 are showed.
#' show(head(jacob$Jacobian.net2, 10))
#' @author Julio E. Sandubete, Lorenzo Escot
#' @seealso \code{\link{embedding}}
#' @importFrom nnet nnet
#' @importFrom NeuralNetTools plotnet
#' @export jacobi
jacobi<-function(x,lag=1,timelapse="FIXED",M0=3,M1=10,H0=2,H1=10,I=100,pre.white=TRUE,doplot=TRUE){

  # Checks
  if (is.null(x)){stop("'x' should be a numeric vector, time serie, data frame or matrix depending on the method selected in 'timelapse'")}
  if (lag < 1){stop("wrong value for the reconstruction delay")}
  if (is.null(timelapse)){stop("'timelapse' should be 'FIXED' or 'VARIABLE'")}
  if ((M0 < 3) | (M1 < M0) | M1 > length(x)){stop("wrong value for the embedding dimension")}
  if ((H0 < 2) | (H1 < H0) | H1 > length(x)){stop("wrong value for the number of neurones in the hidden layer")}
  if (I<1){stop("wrong value of neural networks iterations")}
  if (M0 == M1){
    net.fit<-matrix(nrow = 1,ncol = 3,dimnames = list("Network",c("Embedding","Neurones","Bayesian IC")))
  } else {
    net.fit<-matrix(nrow = M1-(M0-1),ncol = 3,dimnames = list(paste(rep("Network",M1-(M0-1)),c(1:(M1-(M0-1)))),c("Embedding","Neurones","Bayesian IC")))
  }
  if (H0 == H1){
    opt<-matrix(nrow = 1,ncol = 3,dimnames = list("Network",c("Embedding","Neurones","Bayesian IC")))
  } else {
    opt<-matrix(nrow = H1-(H0-1),ncol = 3,dimnames = list(paste(rep("Network",H1-(H0-1)),c(1:(H1-(H0-1)))),c("Embedding","Neurones","Bayesian IC")))
  }

  # Settings
  BIC<-c()
  p<-c() #new
  f<-c() #new
  c<-0
  for (m in M0:M1){
    b<-0
    c<-c+1
    embd<-embedding(x,m=m,lag=lag,timelapse = timelapse)
    y<-as.data.frame(embd[,1])
    xo<-embd[,c(2:m)]
    n<-nrow(xo)
    for (h in H0:H1){
      b<-b+1
      w0<-w0inicial(xo,y,I=I,h=h,m=m,seed.t = T)
      net.nn<-nnet::nnet(xo,y,size=h,maxit=1e5,linout = T,trace=T)
      p[b]<-length(net.nn$wts)
      f[b]<-n*log(net.nn$value/n)
      BIC[b]<-p[b]*log(n)+f[b]
      opt[b,]<-cbind(m,h,BIC[b])
    }
    net.fit[c,]<-opt[which.min(opt[,3]),]
  }
  derivatives<-vector("list",nrow(net.fit))

  # Neural networks with a better fit
  for(k in 1:nrow(net.fit)){
    m<-as.numeric(net.fit[k,1])
    h<-as.numeric(net.fit[k,2])
    embd<-embedding(x,m=m,lag=lag,timelapse = timelapse)
    y<-as.data.frame(embd[,1])
    xo<-embd[,c(2:m)]
    w0<-w0inicial(xo,y,I=I,h=h,m=m,seed.t = T)
    net.nn<-nnet::nnet(xo,y,size=h,Wts=w0, maxit=500,linout = T,trace=F)

    # Plots
    if (doplot==TRUE){
      net.nn$call$x<-xo
      net.nn$call$y<-y
      NeuralNetTools::plotnet(net.nn,y_names="Yt",pos.col="black",neg.col="grey",alpha.val=0.7,circle.cex=4,pad_x=.9)
    }

    # Calculation and evaluation of partial derivatives
    w<-net.nn$wts
    if (pre.white==TRUE) {
      yo.net<-as.vector(net.nn$fitted.values)
      xo.net<-embedding(yo.net,m=m,lag=1,timelapse="FIXED")
      xo.net<-xo.net[,c(2:m)]
    } else {
      xo.net<-xo
    }
    xo.net<-as.matrix(xo.net)
    z<-matrix(nrow = dim(xo.net)[1],ncol=h)
    dphi<-matrix(nrow = dim(xo.net)[1], ncol=h)
    zout<-matrix(nrow = dim(xo.net)[1],ncol=h)
    dzout<-matrix(nrow = dim(xo.net)[1],ncol=(m-1))
    wout<-c()
    wz<-c()
    for (i in 1:h){
      z[,i]<-w[i+(i-1)*(m-1)]+xo.net%*%w[(i+(i-1)*(m-1)+1):((i+(i-1)*(m-1))+(m-1))]
      dphi[,i]<-exp(z[,i])/(1+exp(z[,i]))^2
      wout[i]<-w[1+i+h*m]
      zout[,i]<-dphi[,i]*wout[i]
    }
    for (j in 1:(m-1)){
      for(s in 1:h){
        wz[s]<-w[(s-1)*m+(j+1)]
        wzout<-cbind(wz)
      }
      dzout[,j]<-zout%*%wzout
    }
    dzout<-as.data.frame(dzout)
    name.c<-paste(make.unique(rep("Derivative Xt-", m), sep = ""),"lag",sep="")
    name.r<-make.unique(rep("t=", (nrow(dzout)+1)), sep = "")
    colnames(dzout)[c(1:(m-1))]<-name.c[c(2:m)]
    rownames(dzout)<-name.r[c(2:(nrow(dzout)+1))]
    derivatives[[k]]<-dzout
  }

  # Output
  return(c(list(Network.set=net.fit),Jacobian.net=derivatives))
}

################################################################################
# Internal function to estimate the initial conditions of each neural network
w0inicial<-function(x,y,rangx=1/max(abs(x)),I=100,h=2,m=2,seed.t=TRUE,seed=56666459){

  # Settings
  Wts0<-NULL
  ECMWts0=NULL
  if(seed.t==T){
    set.seed(seed)
  } else {
    set.seed()
  }

  # Estimation of the initial conditions
  for (i in 1:I){
    W0<-runif(h*m+h+1, min=-rangx, max=rangx)
    net.nn<-nnet::nnet(x,y,size=h,Wts=W0,maxit=0,linout = T,trace=F)
    if(i==1){
      ECMWts0<-net.nn$value
      Wts0<-W0
    }
    if(net.nn$value<ECMWts0){
      ECMWts0<-net.nn$value
      Wts0<-W0
    }
  }

  # Output
  return(Wts0)
}





