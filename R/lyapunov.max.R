################################################################################
###### DChaos: an R package for detecting chaotic signals in time series #######
################################################################################
# FUNCTION:         DESCRIPTION:
#   lyapunov.max      Estimates the largest Lyapunov Exponent (NORMA-2)
################################################################################
#' Estimation of the largest Lyapunov Exponent
#' @name lyapunov.max
#' @aliases lyapunov.max
#' @description
#' This function estimates the largest Lyapunov exponent through the Norma-2 method considering the argument set selected by the user.
#' @param x a matrix or data frame containing the partial derivatives of jacobian.
#' @param blocking a character denoting if you consider the full sample \code{FULL}, the non-overlapping sample \code{NOVER}, the equally spaced sample \code{EQS} or the bootstrap sample \code{BOOT} (Default \code{FULL}).
#' @param B a non-negative integer denoting the number of bootstrap iterations (Default 100).
#' @param doplot a logical value denoting if you want to draw a plot \code{TRUE} or not \code{FALSE}. If it is \code{TRUE} shows as many graphs as embedding dimensions have been considered. Each of them represents the evolution of the largest Lyapunov exponent values for the whole period considering the blocking method selected by the user (Default \code{TRUE}).
#' @return A list containing the largest Lyapunov exponent considering the parameter set selected by the user. The dataset saved by each blocking method are the estimated Lyapunov exponent value, the standard error, the z-test value and the p-value. If the user chooses the non-overlapping sample, the equally spaced sample or the bootstrap sample the mean and median values of the Lyapunov exponent are showed. Also some details about the embedding dimension, the sample size, the block length and the block number are recorded.
#' @references Ellner, S., Gallant, A., McCaffrey, D., Nychka, D. 1991 Convergence rates and data requirements for jacobian-based estimates of lyapunov exponents from data. Physics Letters A 153(6):357-363.
#' @references McCaffrey, D.F., Ellner, S., Gallant, A.R., Nychka, D.W. 1992 Estimating the lyapunov exponent of a chaotic system with nonparametric regression. Journal of the American Statistical Association 87(419):682-695.
#' @references Nychka, D., Ellner, S., Gallant, A.R., McCaffrey, D. 1992 Finding chaos in noisy systems. Journal of the Royal Statistical Society 54(2):399-426.
#' @references Whang, Y.J., Linton, O. 1999 The asymptotic distribution of nonparametric estimates of the lyapunov exponent for stochastic time series. Journal of Econometrics 91(1):1-42.
#' @references Shintani, M., Linton, O. 2004 Nonparametric neural network estimation of Lyapunov exponents and a direct test for chaos. Journal of Econometrics 120(1):1-33.
#' ## We show below an example considering time series from the
#' ## logistic equation. We have estimated the largest Lyapunov
#' ## exponent considering the bootstrap sample for an embedding
#' ## dimension (m=4). First of all, we need to estimates the
#' ## partial derivatives of the jacobian.
#' data<-logistic.ts(u.min=4,u.max=4,B=100,doplot=FALSE)
#' ts<-data$`Logistic 100`$time.serie
#' jacob<-jacobi(ts,M0=4,M1=4,doplot=FALSE)
#' deriv<-jacob$Jacobian.net
#' lyapu<-lyapunov.max(deriv,blocking="BOOT",B=100)
#' show(lyapu$Exponent)
#' @author Julio E. Sandubete, Lorenzo Escot
#' @seealso \code{\link{lyapunov.spec}}, \code{\link{lyapunov}}
#' @importFrom pracma normest
#' @importFrom pracma zeros
#' @importFrom sandwich lrvar
#' @importFrom stats pnorm
#' @importFrom graphics abline
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom stats median
#' @export lyapunov.max
lyapunov.max<-function(x,blocking=c("FULL","NOVER","EQS","BOOT"),B=100,doplot=TRUE){

  # Checks
  if (is.null(x)){stop("'x' should be a matrix or data frame containing the partial derivatives of the jacobian")}
  if (is.matrix(x) | is.data.frame(x)){
  if (is.null(blocking)){stop("'blocking' should be 'FULL', 'NOVER', 'EQS' or 'BOOT'")}
  if (B<1){stop("wrong value of bootstrap iterations")}

  # Settings
  N<-nrow(x)
  m<-ncol(x)
  blocking = match.arg(blocking)

  # Estimates the largest Lyapunov Exponent: Full sample
  lyap_max_full<-function(x,doplot=TRUE){

    # Lyapunov exponents
    lvmax<-matrix(0,N,1)
    J<-matrix(0,m,m)
    TM<-diag(1,m)
    v<-matrix(0,m,1)
    v[1,1]<-1
    for (i in 1:N){
      J<-rbind(as.numeric(x[i,1:m]),diag(1,m-1,m))
      if (any(is.nan(J[1,])) | any(is.infinite(J[1,]))){
        lvmax[i]<-lvmax[i-1]
      } else {
        TM<-J%*%TM
        # To avoid that the value of Jacobians tend to infinity after multiplying n-times
        if (any(is.nan(TM[1,])) | any(is.infinite(TM[1,]))){
          lpv<-lvmax[c((i-101):(i-1))]
          k<-rep(seq(1,100,1),ceiling((N-i)/100))
          xo<-i
          s<-1
          for (j in xo:N){
            z<-k[s]
            lvmax[j]<-lpv[z]-runif(1,0.0001,0.001)
            s<-s+1
          }
          break
        } else {
          lvmax[i]<-(1/i)*log(max(svd(TM%*%v)$d))
        }
      }
    }

    # Plots
    if (doplot==TRUE){
      par(mfrow=c(1,1))
      plot(lvmax[c(100:N)],type="l",xlab = "Number of observations", ylab = "Values",main=paste("Largest Lyapunov exponent\n","Full sample | m","=",m,""),col="darkgray")
      abline(h=0,col="steelblue")
    }

    # Eta errors
    etalvmax<-pracma::zeros(N,1)
    etalvmax[1]<-lvmax[1]-lvmax[N]
    for (j in 2:N){
      etalvmax[j]<-j*lvmax[j]-(j-1)*lvmax[j-1]-lvmax[N]
    }

    # Long-run variance of the eta's mean
    varmax<-(sandwich::lrvar(etalvmax, type = c("Andrews"), prewhite = FALSE, adjust = FALSE,kernel = c("Quadratic Spectral"), aprox=c("ARMA(1,1)"),verbose=F))^0.5

    # Hypothesis contrast
    Ztestmax<-lvmax[N]/(varmax/sqrt(N))
    p.value.max<-pnorm(Ztestmax)

    # Output
    LE<-list("Exponent"=matrix(c(lvmax[N],varmax,Ztestmax,p.value.max), nrow=1,ncol=4, byrow=F, dimnames=list("Exponent",c("Estimate","Std. Error","z value","Pr(>|z|)"))),
             "Embedding"=m+1,"Sample"=N
            )
    return(LE)
  }

  if (blocking == "FULL")
    LE = lyap_max_full(x,doplot=doplot)

  # Estimates the largest Lyapunov Exponent: Non-overlapping sample
  lyap_max_over<-function(x,doplot=TRUE){

    # Lyapunov exponents
    M<-trunc(36.2*(N/log(N))^(1/6))
    B<-floor(N/M)
    lvmax<-matrix(0,M,B)
    for (b in 1:B){
      J<-matrix(0,m,m)
      TM<-diag(1,m)
      v<-matrix(0,m,1)
      v[1,1]<-1
      for (i in (1+(b-1)*M):(b*M)){
        J<-rbind(as.numeric(x[i,1:m]),diag(1,m-1,m))
        if (any(is.nan(J[1,])) | any(is.infinite(J[1,]))){
          lvmax[i-(b-1)*M,b]<-lvmax[i-1-(b-1)*M,b]
        } else {
          TM<-J%*%TM
          lvmax[i-(b-1)*M,b]<-(1/(i-(b-1)*M))*log(max(svd(TM%*%v)$d))
        }
      }
    }

    # Plots
    if (doplot==TRUE){
      par(mfrow=c(1,1))
      plot(lvmax[,1],ylim=c(min(lvmax[1:M,]),max(lvmax[1:M,])),type="l",xlab="Block length", ylab = "Values",main=paste("Largest Lyapunov exponent\n","Non-overlapping sample | m","=",m,""),col="darkgray")
      abline(h=0,col="steelblue")
      for (b in 2:B){
        lines(lvmax[,b],type="l",col="darkgray")
      }
    }

    # Eta errors
    etalvmax<-pracma::zeros(M,B)
    for (b in 1:B){
      etalvmax[1,b]<-lvmax[1,b]-lvmax[M,b]
      for (j in 2:M){
        etalvmax[j,b]<-j*lvmax[j,b]-(j-1)*lvmax[j-1,b]-lvmax[M,b]
      }
    }

    # Long-run variance of the eta's mean
    varmax<-pracma::zeros(1,B)
    for (b in 1:B){
      varmax[1,b]<-(sandwich::lrvar(etalvmax[,b], type = c("Andrews"), prewhite = FALSE, adjust = FALSE,kernel = c("Quadratic Spectral"), aprox=c("ARMA(1,1)"),verbose=FALSE))^0.5
    }

    # Hypothesis contrast
    # Mean value
    lvmax_mean<-mean(lvmax[M,])
    sdlvmax_mean<-mean(varmax[1,])
    Ztestmax_mean<-lvmax_mean/(sdlvmax_mean/sqrt(M))
    p.value.max_mean<-pnorm(Ztestmax_mean)

    # Median value
    lvmax_median<-median(lvmax[M,])
    sdlvmax_median<-median(varmax[1,])
    Ztestmax_median<-lvmax_median/(sdlvmax_median/sqrt(M))
    p.value.max_median<-pnorm(Ztestmax_median)

    # Output
    LE<-list("Exponent"=matrix(c(lvmax_mean,lvmax_median,sdlvmax_mean,sdlvmax_median,Ztestmax_mean,Ztestmax_median,p.value.max_mean,p.value.max_median),nrow=2,ncol=4, byrow=F, dimnames=list(c("Exponent-Mean","Exponent-Median"),c("Estimate","Std. Error","z value","Pr(>|z|)"))),
             "Embedding"=m+1,"Sample"=N,"Block.number"=B,"Block.size"=M
             )
    return(LE)
  }

  if (blocking == "NOVER")
    LE = lyap_max_over(x,doplot=doplot)

  # Estimates the largest Lyapunov Exponent: Equally spaced sample
  lyap_max_ess<-function(x,doplot=TRUE){

    # Lyapunov exponents
    M<-trunc(36.2*(N/log(N))^(1/6))
    B<-floor(N/M)
    lvmax<-matrix(0,M,B)
    for (b in 1:B){
      J<-matrix(0,m,m)
      TM<-diag(1,m)
      v<-matrix(0,m,1)
      v[1,1]<-1
      j<-0
      for (i in seq(from=b,to=(B*M), by=B)){
        j<-j+1
        J<-rbind(as.numeric(x[i,1:m]),diag(1,m-1,m))
        if (any(is.nan(J[1,])) | any(is.infinite(J[1,]))){
          lvmax[j,b]<-lvmax[j-1,b]
        } else {
          TM<-J%*%TM
          lvmax[j,b]<-(1/j)*log(max(svd(TM%*%v)$d))
        }
      }
    }

    # Plots
    if (doplot==TRUE){
      par(mfrow=c(1,1))
      plot(lvmax[,1],ylim=c(min(lvmax[1:M,]),max(lvmax[1:M,])),type="l",xlab="Block length", ylab = "Values",main=paste("Largest Lyapunov exponent\n","Equally spaced sample | m","=",m,""),col="darkgray")
      abline(h=0,col="steelblue")
      for (b in 2:B){
        lines(lvmax[,b],type="l",col="darkgray")
      }
    }

    # Eta errors
    etalvmax<-pracma::zeros(M,B)
    for (b in 1:B){
      etalvmax[1,b]<-lvmax[1,b]-lvmax[M,b]
      for (j in 2:M){
        etalvmax[j,b]<-j*lvmax[j,b]-(j-1)*lvmax[j-1,b]-lvmax[M,b]
      }
    }

    # Long-run variance of the eta's mean
    varmax<-pracma::zeros(1,B)
    for (b in 1:B){
      varmax[1,b]<-(sandwich::lrvar(etalvmax[,b], type = c("Andrews"), prewhite = FALSE, adjust = FALSE,kernel = c("Quadratic Spectral"), aprox=c("ARMA(1,1)"),verbose=FALSE))^0.5
    }

    # Hypothesis contrast
    # Mean value
    lvmax_mean<-mean(lvmax[M,])
    sdlvmax_mean<-mean(varmax[1,])
    Ztestmax_mean<-lvmax_mean/(sdlvmax_mean/sqrt(M))
    p.value.max_mean<-pnorm(Ztestmax_mean)

    # Median value
    lvmax_median<-median(lvmax[M,])
    sdlvmax_median<-median(varmax[1,])
    Ztestmax_median<-lvmax_median/(sdlvmax_median/sqrt(M))
    p.value.max_median<-pnorm(Ztestmax_median)

    # Output
    LE<-list("Exponent"=matrix(c(lvmax_mean,lvmax_median,sdlvmax_mean,sdlvmax_median,Ztestmax_mean,Ztestmax_median,p.value.max_mean,p.value.max_median),nrow=2,ncol=4, byrow=F, dimnames=list(c("Exponent-Mean","Exponent-Median"),c("Estimate","Std. Error","z value","Pr(>|z|)"))),
             "Embedding"=m+1,"Sample"=N,"Block.number"=B,"Block.size"=M
            )
    return(LE)
  }

  if (blocking == "EQS")
    LE = lyap_max_ess(x,doplot=doplot)

  # Estimates the largest Lyapunov Exponent: Bootstrap sample
  lyap_max_boot<-function(x,B=100,doplot=TRUE){

    # Lyapunov exponents
    M<-trunc(36.2*(N/log(N))^(1/6))
    B<-B
    lvmax<-matrix(0,M,B)
    for (b in 1:B){
      J<-matrix(0,m,m)
      TM<-diag(1,m)
      v<-matrix(0,m,1)
      v[1,1]<-1
      j<-0
      mboostsample <- sort(sample(N, M,replace=FALSE))
      for (i in mboostsample){
        j<-j+1
        J<-rbind(as.numeric(x[i,1:m]),diag(1,m-1,m))
        if (any(is.nan(J[1,])) | any(is.infinite(J[1,]))){
          lvmax[j,b]<-lvmax[j-1,b]
        } else {
          TM<-J%*%TM
          lvmax[j,b]<-(1/j)*log(max(svd(TM%*%v)$d))
        }
      }
    }

    # Plots
    if (doplot==TRUE){
      par(mfrow=c(1,1))
      lvmax.mean<-c()
      plot(lvmax[,1],ylim=c(min(lvmax[1:M,]),max(lvmax[1:M,])),type="l",xlab="Block length", ylab = "Values",main=paste("Largest Lyapunov exponent\n","Bootstrap sample | m","=",m,""),col="darkgray")
      abline(h=0,col="steelblue")
      for (b in 2:B){
        lines(lvmax[,b],type="l",col="darkgray")
      }
      for (a in 1:M){
        lvmax.mean[a]<-mean(lvmax[a,c(1:B)])
      }
      lines(lvmax.mean,type="l",col="indianred1")
    }

    # Eta errors
    etalvmax<-pracma::zeros(M,B)
    for (b in 1:B){
      etalvmax[1,b]<-lvmax[1,b]-lvmax[M,b]
      for (j in 2:M){
        etalvmax[j,b]<-j*lvmax[j,b]-(j-1)*lvmax[j-1,b]-lvmax[M,b]
      }
    }

    # Long-run variance of the eta's mean
    varmax<-pracma::zeros(1,B)
    for (b in 1:B){
      varmax[1,b]<-(sandwich::lrvar(etalvmax[,b], type = c("Andrews"), prewhite = FALSE, adjust = FALSE,kernel = c("Quadratic Spectral"), aprox=c("ARMA(1,1)"),verbose=FALSE))^0.5
    }

    # Hypothesis contrast
    # Mean value
    lvmax_mean<-mean(lvmax[M,abs(outliers::scores(lvmax[M,], type = c("iqr")))<1.5])
    sdlvmax_mean<-mean(varmax[1,abs(outliers::scores(varmax[1,], type = c("iqr")))<1.5])
    Ztestmax_mean<-lvmax_mean/(sdlvmax_mean/sqrt(M))
    p.value.max_mean<-pnorm(Ztestmax_mean)

    # Median value
    lvmax_median<-median(lvmax[M,])
    sdlvmax_median<-median(varmax[1,])
    Ztestmax_median<-lvmax_median/(sdlvmax_median/sqrt(M))
    p.value.max_median<-pnorm(Ztestmax_median)

    # Output
    LE<-list("Exponent"=matrix(c(lvmax_mean,lvmax_median,sdlvmax_mean,sdlvmax_median,Ztestmax_mean,Ztestmax_median,p.value.max_mean,p.value.max_median),nrow=2,ncol=4, byrow=F, dimnames=list(c("Exponent-Mean","Exponent-Median"),c("Estimate","Std. Error","z value","Pr(>|z|)"))),
             "Embedding"=m+1,"Sample"=N,"Block.number"=B,"Block.size"=M
            )
    return(LE)
  }

  if (blocking == "BOOT")
    LE = lyap_max_boot(x,B=B,doplot=doplot)

  # Output
  return(LE)
  } else {
  stop("'x' should be a matrix or data frame containing the partial derivatives of the jacobian")
  }
}



