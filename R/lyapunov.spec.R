################################################################################
###### DChaos: an R package for detecting chaotic signals in time series #######
################################################################################
# FUNCTION:         DESCRIPTION:
#   lyapunov.spec     Estimates the Lyapunov exponent spectrum (QR decomposition)
################################################################################
#' Estimation of the Lyapunov exponent spectrum
#' @name lyapunov.spec
#' @aliases lyapunov.spec
#' @description
#' This function estimates the Lyapunov exponent spectrum through the QR decomposition method considering the argument set selected by the user.
#' @param x a matrix or data frame containing the partial derivatives of jacobian.
#' @param blocking a character denoting if you consider the full sample \code{FULL}, the non-overlapping sample \code{NOVER}, the equally spaced sample \code{EQS} or the bootstrap sample \code{BOOT} (Default \code{FULL}).
#' @param B a non-negative integer denoting the number of bootstrap iterations (Default 100).
#' @param doplot a logical value denoting if you want to draw a plot \code{TRUE} or not \code{FALSE}. If it is \code{TRUE} shows as many graphs as embedding dimensions have been considered. Each of them represents the evolution of the Lyapunov exponent spectrum values for the whole period considering the blocking method selected by the user (Default \code{TRUE}).
#' @return A list containing the Lyapunov exponent spectrum considering the parameter set selected by the user. The dataset saved by each blocking method are the estimated Lyapunov exponent value, the standard error, the z-test value and the p-value. If the user chooses the non-overlapping sample, the equally spaced sample or the bootstrap sample the mean and median values of the Lyapunov exponent are showed. Also some details about the embedding dimension, the sample size, the block length and the block number are recorded.
#' @references Ellner, S., Gallant, A., McCaffrey, D., Nychka, D. 1991 Convergence rates and data requirements for jacobian-based estimates of lyapunov exponents from data. Physics Letters A 153(6):357-363.
#' @references McCaffrey, D.F., Ellner, S., Gallant, A.R., Nychka, D.W. 1992 Estimating the lyapunov exponent of a chaotic system with nonparametric regression. Journal of the American Statistical Association 87(419):682-695.
#' @references Nychka, D., Ellner, S., Gallant, A.R., McCaffrey, D. 1992 Finding chaos in noisy systems. Journal of the Royal Statistical Society 54(2):399-426.
#' @references Whang, Y.J., Linton, O. 1999 The asymptotic distribution of nonparametric estimates of the lyapunov exponent for stochastic time series. Journal of Econometrics 91(1):1-42.
#' @references Shintani, M., Linton, O. 2004 Nonparametric neural network estimation of Lyapunov exponents and a direct test for chaos. Journal of Econometrics 120(1):1-33.
#' @examples
#' ## We show below an example considering time series from the
#' ## logistic equation. We have estimated the Lyapunov exponent
#' ## spectrum considering the bootstrap sample for an embedding
#' ## dimension (m=4). First of all, we need to estimates the
#' ## partial derivatives of the jacobian.
#' data<-logistic.ts(u.min=4,u.max=4,B=100,doplot=FALSE)
#' ts<-data$`Logistic 100`$time.serie
#' jacob<-jacobi(ts,M0=4,M1=4,doplot=FALSE)
#' deriv<-jacob$Jacobian.net
#' lyapu<-lyapunov.spec(deriv,blocking="BOOT",B=10)
#' show(lyapu$Exponent)
#' @author Julio E. Sandubete, Lorenzo Escot
#' @seealso \code{\link{lyapunov.max}}, \code{\link{lyapunov}}
#' @importFrom pracma normest
#' @importFrom pracma zeros
#' @importFrom sandwich lrvar
#' @importFrom stats pnorm
#' @importFrom graphics abline
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom stats median
#' @export lyapunov.spec
lyapunov.spec<-function(x,blocking=c("FULL","NOVER","EQS","BOOT"),B=100,doplot=TRUE){

  # Checks
  if (is.null(x)){stop("'x' should be a matrix or data frame containing the partial derivatives of the jacobian")}
  if (is.matrix(x) | is.data.frame(x)){
  if (is.null(blocking)){stop("'blocking' should be 'FULL', 'NOVER', 'EQS' or 'BOOT'")}
  if (B<1){stop("wrong value of bootstrap iterations")}

  # Settings
  N<-nrow(x)
  m<-ncol(x)
  blocking = match.arg(blocking)

  # Estimates the Lyapunov exponent spectrum: Full sample
  lyap_spec_full<-function(x,doplot=TRUE){

    # Lyapunov exponents
    lv<-matrix(0,N,m)
    lpv<-matrix(0,m,1)
    J<-matrix(0,m,m)
    TM<-diag(1,m)
    for (i in 1:N){
      J<-rbind(as.numeric(x[i,1:m]),diag(1,m-1,m))
      TM<-J%*%TM
      qro<-qr(TM)
      Q<-qr.Q(qro)
      R<-qr.R(qro)
      Q<-Q%*%solve(diag(diag(sign(R))))
      R<-R%*%diag(diag(sign(R)))
      for (j in 1:m){
        lpv[j]<-lpv[j]+log(R[j,j])
        lv[i,j]<-lpv[j]/i
      }
      TM<-Q
    }
    lpv<-lpv/N

    # Plots
    if (doplot==TRUE){
      if (m < 5){
        par(mfrow=c(ceiling(m/2),2))
      } else {
        par(mfrow=c(ceiling(m/3),3))
      }
      for (r in 1:m){
        plot(lv[,r], type="l",xlab="Number of observations",ylab="Values",main=paste("Lyapunov exponent",r,"\n","Full sample | m","=",m,""),col="darkgray")
        abline(h=0,col="steelblue")
      }
      par(mfrow=c(1,1))
    }

    # Eta errors
    etalv<-pracma::zeros(N,m)
    for (r in 1:m){
      etalv[1,r]<-lv[1,r]-lv[N,r]
    }
    for (j in 2:N){
      for (r in 1:m){
        etalv[j,r]<-j*lv[j,r]-(j-1)*lv[j-1,r]-lv[N,r]
      }
    }

    # Long-run variance of the eta's mean
    varlv2<-pracma::zeros(1,m)
    varlv2<-diag((sandwich::lrvar(etalv, type = c("Andrews"), prewhite = FALSE, adjust = FALSE,kernel = c("Quadratic Spectral"), aprox=c("ARMA(1,1)"),verbose=F)))^0.5

    # Hypothesis contrast (H0: λ >= 0)
    Ztest<-pracma::zeros(1,m)
    p.value<-pracma::zeros(1,m)
    for (r in 1:m){
      Ztest[r]<-lpv[r]/(varlv2[r]/sqrt(N))
      p.value[r]<-pnorm(Ztest[r])
    }

    # Output
    LE<-list("Exponent"=matrix(c(lv[N,],varlv2,Ztest,p.value), nrow=m,ncol=4, byrow=F, dimnames=list(paste("Exponent", seq(1:m)),c("Estimate","Std. Error","z value","Pr(>|z|)"))),
             "Embedding"=m+1,"Sample"=N
            )
    return(LE)
  }

  if (blocking == "FULL")
    LE = lyap_spec_full(x,doplot=doplot)

  # Estimates the Lyapunov exponent spectrum: Non-overlapping sample
  lyap_spec_over<-function(x,doplot=TRUE){

    # Lyapunov exponents
    M<-trunc(36.2*(N/log(N))^(1/6))
    B<-floor(N/M)
    lpv<-array(0,dim=c(B,m))
    lv<-array(0,dim=c(M,B,m))
    for (b in 1:B){
      J<-matrix(0,m,m)
      TM<-diag(1,m)
      for (i in (1+(b-1)*M):(b*M)){
        J<-rbind(as.numeric(x[i,1:m]),diag(1,m-1,m))
        TM<-J%*%TM
        qro<-qr(TM)
        Q<-qr.Q(qro)
        R<-qr.R(qro)
        Q<-Q%*%solve(diag(diag(sign(R))))
        R<-R%*%diag(diag(sign(R)))
        for (j in 1:m){
          lpv[b,j]<-lpv[b,j]+log(R[j,j])
          lv[i-(b-1)*M,b,j]<-lpv[b,j]/(i-(b-1)*M)
        }
        TM<-Q
      }
    }
    lpv<-lpv/M

    # Plots
    if (doplot==TRUE){
      if (m < 5){
        par(mfrow=c(ceiling(m/2),2))
      } else {
        par(mfrow=c(ceiling(m/3),3))
      }
      for (r in 1:m){
        plot(lv[,1,r],ylim=c(min(lv[,,r]),max(lv[,,r])), type="l",xlab="Block length",ylab="Values", main=paste("Lyapunov exponent",r,"\n","Non-overlapping sample | m","=",m,""),col="darkgray")
        abline(h=0,col="steelblue")
        for (b in 2:B){
          lines(lv[,b,r], type="l",col="darkgray")
        }
      }
      par(mfrow=c(1,1))
    }

    # Eta errors
    etalv<-array(0,dim=c(M,B,m))
    for (b in 1:B){
      for (r in 1:m){
        etalv[1,b,r]<-lv[1,b,r]-lv[M,b,r]
      }
      for (j in 2:M){
        for (r in 1:m){
          etalv[j,b,r]<-j*lv[j,b,r]-(j-1)*lv[j-1,b,r]-lv[M,b,r]
        }
      }
    }

    # Long-run variance of the eta's mean
    varlv2<-array(0,dim=c(B,m))
    for (b in 1:B){
      varlv2[b,]<-diag((sandwich::lrvar(etalv[,b,], type = c("Andrews"), prewhite = FALSE, adjust = FALSE,kernel = c("Quadratic Spectral"), aprox=c("ARMA(1,1)"),verbose=FALSE)))^0.5
    }

    # Hypothesis contrast (H0: λ >= 0)
    # Mean value
    lpv_mean<-apply(lpv,2,mean)
    sdlpv_mean<-apply(varlv2,2,mean)
    Ztest_mean<-pracma::zeros(1,m)
    p.value_mean<-pracma::zeros(1,m)
    for (r in 1:m){
      Ztest_mean[r]<-lpv_mean[r]/(sdlpv_mean[r]/sqrt(M))
      p.value_mean[r]<-pnorm(Ztest_mean[r])
    }

    # Median value
    lpv_median<-apply(lpv,2,median)
    sdlpv_median<-apply(varlv2,2,median)
    Ztest_median<-pracma::zeros(1,m)
    p.value_median<-pracma::zeros(1,m)
    for (r in 1:m){
      Ztest_median[r]<-lpv_median[r]/(sdlpv_median[r]/sqrt(M))
      p.value_median[r]<-pnorm(Ztest_median[r])
    }

    # Output
    LE<-list("Exponent"=rbind(matrix(c(lpv_mean,sdlpv_mean,Ztest_mean,p.value_mean), nrow=m,ncol=4, byrow=F, dimnames=list(paste("Exponent-Mean", seq(1:m)),c("Estimate","Std. Error","z value","Pr(>|z|)"))),
                              matrix(c(lpv_median,sdlpv_median,Ztest_median,p.value_median), nrow=m,ncol=4, byrow=F, dimnames=list(paste("Exponent-Median", seq(1:m)),c("Estimate","Std. Error","z value","Pr(>|z|)")))),
             "Embedding"=m+1,"Sample"=N,"Blocks.number"=B,"Block.size"=M
            )
    return(LE)
  }

  if (blocking == "NOVER")
    LE = lyap_spec_over(x,doplot=doplot)

  # Estimates the Lyapunov exponent spectrum: Equally spaced sample
  lyap_spec_ess<-function(x,doplot=TRUE){

    # Lyapunov exponents
    M<-trunc(36.2*(N/log(N))^(1/6))
    B<-floor(N/M)
    lpv<-array(0,dim=c(B,m))
    lv<-array(0,dim=c(M,B,m))
    for (b in 1:B){
      J<-matrix(0,m,m)
      TM<-diag(1,m)
      j<-0
      for (i in seq(from=b,to=(B*M), by=B)){
        j<-j+1
        J<-rbind(as.numeric(x[i,1:m]),diag(1,m-1,m))
        TM<-J%*%TM
        qro<-qr(TM)
        Q<-qr.Q(qro)
        R<-qr.R(qro)
        Q<-Q%*%solve(diag(diag(sign(R))))
        R<-R%*%diag(diag(sign(R)))
        for (r in 1:m){
          lpv[b,r]<-lpv[b,r]+log(R[r,r])
          lv[j,b,r]<-lpv[b,r]/(j)
        }
        TM<-Q
      }
    }
    lpv<-lpv/M

    # Plots
    if (doplot==TRUE){
      if (m < 5){
        par(mfrow=c(ceiling(m/2),2))
      } else {
        par(mfrow=c(ceiling(m/3),3))
      }
      for (r in 1:m){
        plot(lv[,1,r], ylim=c(min(lv[,,r]),max(lv[,,r])), type="l",xlab="Block length", ylab="Values", main=paste("Lyapunov exponent",r,"\n","Equally spaced sample | m","=",m,""),col="darkgray")
        abline(h=0,col="steelblue")
        for (b in 2:B){
          lines(lv[,b,r], type="l",col="darkgray")
        }
      }
      par(mfrow=c(1,1))
    }

    # Eta errors
    etalv<-array(0,dim=c(M,B,m))
    for (b in 1:B){
      for (r in 1:m){
        etalv[1,b,r]<-lv[1,b,r]-lv[M,b,r]
      }
      for (j in 2:M){
        for (r in 1:m){
          etalv[j,b,r]<-j*lv[j,b,r]-(j-1)*lv[j-1,b,r]-lv[M,b,r]
        }
      }
    }

    # Long-run variance of the eta's mean
    varlv2<-array(0,dim=c(B,m))
    for (b in 1:B){
      varlv2[b,]<-diag((sandwich::lrvar(etalv[,b,], type = c("Andrews"), prewhite = FALSE, adjust = FALSE,kernel = c("Quadratic Spectral"), aprox=c("ARMA(1,1)"),verbose=FALSE)))^0.5
    }

    # Hypothesis contrast (H0: λ >= 0)
    # Mean value
    lpv_mean<-apply(lpv,2,mean)
    sdlpv_mean<-apply(varlv2,2,mean)
    Ztest_mean<-pracma::zeros(1,m)
    p.value_mean<-pracma::zeros(1,m)
    for (r in 1:m){
      Ztest_mean[r]<-lpv_mean[r]/(sdlpv_mean[r]/sqrt(M))
      p.value_mean[r]<-pnorm(Ztest_mean[r])
    }

    # Median value
    lpv_median<-apply(lpv,2,median)
    sdlpv_median<-apply(varlv2,2,median)
    Ztest_median<-pracma::zeros(1,m)
    p.value_median<-pracma::zeros(1,m)
    for (r in 1:m){
      Ztest_median[r]<-lpv_median[r]/(sdlpv_median[r]/sqrt(M))
      p.value_median[r]<-pnorm(Ztest_median[r])
    }

    # Output
    LE<-list("Exponent"=rbind(matrix(c(lpv_mean,sdlpv_mean,Ztest_mean,p.value_mean), nrow=m,ncol=4, byrow=F, dimnames=list(paste("Exponent-Mean", seq(1:m)),c("Estimate","Std. Error","z value","Pr(>|z|)"))),
                              matrix(c(lpv_median,sdlpv_median,Ztest_median,p.value_median), nrow=m,ncol=4, byrow=F, dimnames=list(paste("Exponent-Median", seq(1:m)),c("Estimate","Std. Error","z value","Pr(>|z|)")))),
             "Embedding"=m+1,"Sample"=N,"Blocks.number"=B,"Block.size"=M
            )
    return(LE)
  }

  if (blocking == "EQS")
    LE = lyap_spec_ess(x,doplot=doplot)

  # Estimates the Lyapunov exponent spectrum: Bootstrap sample
  lyap_spec_boot<-function(x,B=100,doplot=TRUE){

    # Lyapunov exponents
    M<-trunc(36.2*(N/log(N))^(1/6))
    B<-B
    lpv<-array(0,dim=c(B,m))
    lv<-array(0,dim=c(M,B,m))
    for (b in 1:B){
      J<-matrix(0,m,m)
      TM<-diag(1,m)
      j<-0
      mbootsample <- sort(sample(N, M,replace=FALSE))
      for (i in mbootsample){
        j<-j+1
        J<-rbind(as.numeric(x[i,1:m]),diag(1,m-1,m))
        TM<-J%*%TM
        qro<-qr(TM)
        Q<-qr.Q(qro)
        R<-qr.R(qro)
        Q<-Q%*%solve(diag(diag(sign(R))))
        R<-R%*%diag(diag(sign(R)))
        for (r in 1:m){
          lpv[b,r]<-lpv[b,r]+log(R[r,r])
          lv[j,b,r]<-lpv[b,r]/(j)
        }
        TM<-Q
      }
    }
    lpv<-lpv/M

    # Plots
    if (doplot==TRUE){
      lv.mean<-c()
      if (m < 5){
        par(mfrow=c(ceiling(m/2),2))
      } else {
        par(mfrow=c(ceiling(m/3),3))
      }
      for (r in 1:m){
        plot(lv[,1,r], ylim=c(min(lv[,,r]),max(lv[,,r])), type="l",xlab="Block length", ylab="Values", main=paste("Lyapunov exponent",r,"\n","Bootstrap sample | m","=",m,""),col="darkgray")
        abline(h=0,col="steelblue")
        for (b in 2:B){
          lines(lv[,b,r], type="l",col="darkgray")
        }
        for (a in 1:M){
          lv.mean[a]<-mean(lv[a,c(1:B),r])
        }
        lines(lv.mean,type="l",col="indianred1")
      }
      par(mfrow=c(1,1))
    }

    # Eta errors
    etalv<-array(0,dim=c(M,B,m))
    for (b in 1:B){
      for (r in 1:m){
        etalv[1,b,r]<-lv[1,b,r]-lv[M,b,r]
      }
      for (j in 2:M){
        for (r in 1:m){
          etalv[j,b,r]<-j*lv[j,b,r]-(j-1)*lv[j-1,b,r]-lv[M,b,r]
        }
      }
    }

    # Long-run variance of the eta's mean
    varlv2<-array(0,dim=c(B,m))
    for (b in 1:B){
      varlv2[b,]<-diag((sandwich::lrvar(etalv[,b,], type = c("Andrews"), prewhite = FALSE, adjust = FALSE,kernel = c("Quadratic Spectral"), aprox=c("ARMA(1,1)"),verbose=FALSE)))^0.5
    }

    # Hypothesis contrast (H0: λ >= 0)
    # Mean value
    lpv_mean<-apply(lpv,2,mean)
    sdlpv_mean<-apply(varlv2,2,mean)
    Ztest_mean<-pracma::zeros(1,m)
    p.value_mean<-pracma::zeros(1,m)
    for (r in 1:m){
      Ztest_mean[r]<-lpv_mean[r]/(sdlpv_mean[r]/sqrt(M))
      p.value_mean[r]<-pnorm(Ztest_mean[r])
    }

    # Median value
    lpv_median<-apply(lpv,2,median)
    sdlpv_median<-apply(varlv2,2,median)
    Ztest_median<-pracma::zeros(1,m)
    p.value_median<-pracma::zeros(1,m)
    for (r in 1:m){
      Ztest_median[r]<-lpv_median[r]/(sdlpv_median[r]/sqrt(M))
      p.value_median[r]<-pnorm(Ztest_median[r])
    }

    # Output
    LE<-list("Exponent"=rbind(matrix(c(lpv_mean,sdlpv_mean,Ztest_mean,p.value_mean), nrow=m,ncol=4, byrow=F, dimnames=list(paste("Exponent-Mean", seq(1:m)),c("Estimate","Std. Error","z value","Pr(>|z|)"))),
                              matrix(c(lpv_median,sdlpv_median,Ztest_median,p.value_median), nrow=m,ncol=4, byrow=F, dimnames=list(paste("Exponent-Median", seq(1:m)),c("Estimate","Std. Error","z value","Pr(>|z|)")))),
             "Embedding"=m+1,"Sample"=N,"Blocks.number"=B,"Block.size"=M
            )
    return(LE)
  }

  if (blocking == "BOOT")
    LE = lyap_spec_boot(x,B=B,doplot=doplot)

  # Output
  return(LE)
  } else {
    stop("'x' should be a matrix or data frame containing the partial derivatives of the jacobian")
  }
}



