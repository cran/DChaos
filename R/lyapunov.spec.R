################################################################################
############ DChaos: An R Package for Chaotic Time Series Analysis #############
################################################################################
# FUNCTION:         DESCRIPTION:
#   lyapunov.spec     Estimates the Lyapunov exponent spectrum (QR decomposition)
################################################################################
#' Estimates the Lyapunov exponent spectrum
#' @name lyapunov.spec
#' @aliases lyapunov.spec
#' @description
#' This function estimates the Lyapunov exponent spectrum through the QR decomposition procedure based on the partial derivatives computed by the \code{jacobian.net} function.
#' @param x a \code{matrix}, \code{data.frame} or \code{data.table} containing the partial derivatives computed by the \code{jacobian.net} function.
#' @param blocking a character denoting the blocking method chosen for figuring out the Lyapunov exponent spectrum through the QR decomposition procedure. Available options are \code{FULL} if the user considers the full sample, \code{NOVER} if the user considers the non-overlapping sample, \code{EQS} if the user considers the equally spaced sample or \code{BOOT} if the user considers the bootstrap sample (Default \code{BOOT}).
#' @param B a non-negative integer denoting the number of bootstrap iterations (Default 1000).
#' @param doplot a logical value denoting if the user wants to draw a plot \code{TRUE} or not \code{FALSE}. If it is \code{TRUE} the evolution of the Lyapunov exponent values are represented for the whole period considering the blocking method chosen by the user. It shows as many graphs as embedding dimensions have been considered (Default \code{TRUE}).
#' @return This function returns several objects considering the parameter set selected by the user. The Lyapunov exponent spectrum considering the QR decomposition procedure by each blocking method are estimated. It also contains some useful information about the estimated jacobian, the best-fitted feed-forward single hidden layer neural net model, the best set of weights found, the fitted values, the residuals obtained, the best embedding parameters set chosen, the sample size or the block length considered by each blocking method. This function provides the standard error, the z test value and the p-value for testing the null hypothesis \eqn{H0: \lambda_k > 0 for k = 1,2,3, \ldots, m} (full spectrum). Reject the null hypothesis ${H_0}$ means lack of chaotic behaviour. That is, the data-generating process does not have a chaotic attractor because of it does not show the property of sensitivity to initial conditions.
#' @note The DChaos package provides several ways to figure out robustly the neural net estimator of the k-th Lyapunov exponent. On the one hand if the R users have previously obtained the partial derivatives from the \code{jacobian.net} function they can apply directly the function \code{lyapunov.spec} which estimates the Lyapunov exponent spectrum taking into account the QR decomposition procedure. They can also use the function \code{lyapunov.max} which estimates only the largest Lyapunov exponent considering the Norma-2 procedure. Hence the DChaos package allows the R users to choose between two different procedures to obtain the neural net estimator of the k-th Lyapunov exponent and four ways of subsampling by blocks: full sample, non-overlapping sample, equally spaced sample and bootstrap sample. The blocking methods what they do is to split the time-series data into several blocks by estimating a Lyapunov exponent for each subsample. If the R users choose the non-overlapping sample (\code{blocking = "NOVER"}), the equally spaced sample (\code{blocking = "EQS"}) or the bootstrap sample (\code{blocking = "BOOT"}) the mean and median values of the Lyapunov exponent for each block or subsample are saved. By default we recommend using the median value as it is more robust to the presence of outliers. Notice that the parameter \code{B} will only be considered if the R users choose the bootstrap blocking method.
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
#' ## bootstrap blocking method from the best-fitted neural net model and the partial
#' ## derivatives showed in the jacobian.net example.
#' ## jacobian <- DChaos::jacobian.net(data=ts, m=3:3, lag=1:1, timelapse="FIXED", h=2:10)
#' ## summary(jacobian)
#' ## exponent <- DChaos::lyapunov.spec(x=jacobian, blocking="BOOT", B=100, doplot=FALSE)
#' ## summary(exponent)
#' @author Julio E. Sandubete, Lorenzo Escot
#' @importFrom pracma normest
#' @importFrom pracma zeros
#' @importFrom sandwich lrvar
#' @importFrom stats pnorm
#' @importFrom graphics abline
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom stats median
#' @export lyapunov.spec
lyapunov.spec    <- function(x, blocking=c("BOOT","NOVER","EQS","FULL"), B=1000, doplot=TRUE){

  # Checks
  if (is.null(x$jacobian)){stop("'x' should be a matrix or data frame containing the partial derivatives of the jacobian")}
  if (is.matrix(x$jacobian) | is.data.frame(x$jacobian)){
    if (is.null(blocking)){stop("'blocking' should be 'BOOT', 'NOVER', 'EQS' or 'FULL'")}
    if (B<1){stop("wrong value of bootstrap iterations")}

    # Settings
    x <- x$jacobian
    N <- nrow(x)
    m <- ncol(x)
    blocking = match.arg(blocking)

    # Estimates the Lyapunov exponent: Bootstrap sample
    if (blocking == "BOOT"){
      lyap_spec_boot<-function(x,B=1000,doplot=TRUE){

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
            if (m==1){
              if (length(which(diag(sign(R))==0))!=0){
                TM<-Q
                next
              }
              Q<-Q%*%solve(diag(sign(R)))
              R<-R%*%diag(sign(R))
            } else {
              if (length(which(diag(sign(R))==0))!=0){
                TM<-Q
                next
              }
              Q<-Q%*%solve(diag(diag(sign(R))))
              R<-R%*%diag(diag(sign(R)))
            }
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
          if(m==1){
            varlv2[b,]<-(sandwich::lrvar(etalv[,b,], type = c("Andrews"), prewhite = FALSE, adjust = FALSE,kernel = c("Quadratic Spectral"), aprox=c("ARMA(1,1)"),verbose=F))^0.5
          } else {
            varlv2[b,]<-diag((sandwich::lrvar(etalv[,b,], type = c("Andrews"), prewhite = FALSE, adjust = FALSE,kernel = c("Quadratic Spectral"), aprox=c("ARMA(1,1)"),verbose=F)))^0.5
          }
        }

        # Hypothesis contrast
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
        LE<-list("estimator"=c("Lyapunov exponent spectrum"),"procedure"=c("QR decomposition by bootstrap blocking method"),
                 "exponent.mean"=matrix(c(lpv_mean,sdlpv_mean,Ztest_mean,p.value_mean), nrow=m,ncol=4, byrow=F, dimnames=list(paste("Exponent", seq(1:m)),c("Estimate","Std. Error","z value","Pr(>|z|)"))),
                 "exponent.median"=matrix(c(lpv_median,sdlpv_median,Ztest_median,p.value_median), nrow=m,ncol=4, byrow=F, dimnames=list(paste("Exponent", seq(1:m)),c("Estimate","Std. Error","z value","Pr(>|z|)"))),
                 "sample"=N,"block.length"=M,"no.block"=B
        )
        return(LE)
      }
      LE = lyap_spec_boot(x,B=B,doplot=doplot)
    }

    # Estimates the Lyapunov exponent: Non-overlapping sample
    if (blocking == "NOVER"){
      lyap_spec_over<-function(x,doplot=TRUE){

        N<-nrow(x)
        m<-ncol(x)

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
            if (m==1){
              if (length(which(diag(sign(R))==0))!=0){
                TM<-Q
                next
              }
              Q<-Q%*%solve(diag(sign(R)))
              R<-R%*%diag(sign(R))
            } else {
              if (length(which(diag(sign(R))==0))!=0){
                TM<-Q
                next
              }
              Q<-Q%*%solve(diag(diag(sign(R))))
              R<-R%*%diag(diag(sign(R)))
            }
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
          if(m==1){
            varlv2[b,]<-(sandwich::lrvar(etalv[,b,], type = c("Andrews"), prewhite = FALSE, adjust = FALSE,kernel = c("Quadratic Spectral"), aprox=c("ARMA(1,1)"),verbose=F))^0.5
          } else {
            varlv2[b,]<-diag((sandwich::lrvar(etalv[,b,], type = c("Andrews"), prewhite = FALSE, adjust = FALSE,kernel = c("Quadratic Spectral"), aprox=c("ARMA(1,1)"),verbose=F)))^0.5
          }
        }

        # Hypothesis contrast
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
        LE<-list("estimator"=c("Lyapunov exponent spectrum"),"procedure"=c("QR decomposition by non-overlapping blocking method"),
                 "exponent.mean"=matrix(c(lpv_mean,sdlpv_mean,Ztest_mean,p.value_mean), nrow=m,ncol=4, byrow=F, dimnames=list(paste("Exponent", seq(1:m)),c("Estimate","Std. Error","z value","Pr(>|z|)"))),
                 "exponent.median"=matrix(c(lpv_median,sdlpv_median,Ztest_median,p.value_median), nrow=m,ncol=4, byrow=F, dimnames=list(paste("Exponent", seq(1:m)),c("Estimate","Std. Error","z value","Pr(>|z|)"))),
                 "sample"=N,"block.length"=M,"no.block"=B
        )
        return(LE)
      }
      LE = lyap_spec_over(x,doplot=doplot)
    }

    # Estimates the Lyapunov exponent: Equally spaced sample
    if (blocking == "EQS"){
      lyap_spec_ess<-function(x,doplot=TRUE){

        N<-nrow(x)
        m<-ncol(x)

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
            if (m==1){
              if (length(which(diag(sign(R))==0))!=0){
                TM<-Q
                next
              }
              Q<-Q%*%solve(diag(sign(R)))
              R<-R%*%diag(sign(R))
            } else {
              if (length(which(diag(sign(R))==0))!=0){
                TM<-Q
                next
              }
              Q<-Q%*%solve(diag(diag(sign(R))))
              R<-R%*%diag(diag(sign(R)))
            }
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
          if(m==1){
            varlv2[b,]<-(sandwich::lrvar(etalv[,b,], type = c("Andrews"), prewhite = FALSE, adjust = FALSE,kernel = c("Quadratic Spectral"), aprox=c("ARMA(1,1)"),verbose=F))^0.5
          } else {
            varlv2[b,]<-diag((sandwich::lrvar(etalv[,b,], type = c("Andrews"), prewhite = FALSE, adjust = FALSE,kernel = c("Quadratic Spectral"), aprox=c("ARMA(1,1)"),verbose=F)))^0.5
          }
        }

        # Hypothesis contrast
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
        LE<-list("estimator"=c("Lyapunov exponent spectrum"),"procedure"=c("QR decomposition by equally spaced blocking method"),
                 "exponent.mean"=matrix(c(lpv_mean,sdlpv_mean,Ztest_mean,p.value_mean), nrow=m,ncol=4, byrow=F, dimnames=list(paste("Exponent", seq(1:m)),c("Estimate","Std. Error","z value","Pr(>|z|)"))),
                 "exponent.median"=matrix(c(lpv_median,sdlpv_median,Ztest_median,p.value_median), nrow=m,ncol=4, byrow=F, dimnames=list(paste("Exponent", seq(1:m)),c("Estimate","Std. Error","z value","Pr(>|z|)"))),
                 "sample"=N,"block.length"=M,"no.block"=B
        )
        return(LE)
      }
      LE = lyap_spec_ess(x,doplot=doplot)
    }

    # Estimates the Lyapunov exponent: Full sample
    if (blocking == "FULL"){
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
          if (m==1){
            if (length(which(diag(sign(R))==0))!=0){
              TM<-Q
              next
            }
            Q<-Q%*%solve(diag(sign(R)))
            R<-R%*%diag(sign(R))
          } else {
            if (length(which(diag(sign(R))==0))!=0){
              TM<-Q
              next
            }
            Q<-Q%*%solve(diag(diag(sign(R))))
            R<-R%*%diag(diag(sign(R)))
          }
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
        if(m==1){
          varlv2<-(sandwich::lrvar(etalv, type = c("Andrews"), prewhite = FALSE, adjust = FALSE,kernel = c("Quadratic Spectral"), aprox=c("ARMA(1,1)"),verbose=F))^0.5
        } else {
          varlv2<-diag((sandwich::lrvar(etalv, type = c("Andrews"), prewhite = FALSE, adjust = FALSE,kernel = c("Quadratic Spectral"), aprox=c("ARMA(1,1)"),verbose=F)))^0.5
        }

        # Hypothesis contrast
        Ztest<-pracma::zeros(1,m)
        p.value<-pracma::zeros(1,m)
        for (r in 1:m){
          Ztest[r]<-lpv[r]/(varlv2[r]/sqrt(N))
          p.value[r]<-pnorm(Ztest[r])
        }

        # Output
        LE<-list("estimator"=c("Lyapunov exponent spectrum"),"procedure"=c("QR decomposition by full sample method"),
                 "exponent"=matrix(c(lv[N,],varlv2,Ztest,p.value), nrow=m,ncol=4, byrow=F, dimnames=list(paste("Exponent", seq(1:m)),c("Estimate","Std. Error","z value","Pr(>|z|)"))),
                 "sample"=N,"block.length"=N,"no.block"=1
        )
        return(LE)
      }
      LE = lyap_spec_full(x,doplot=doplot)
    }

    # Class definition
    LE <- c(x,LE,nprint=0)
    class(LE) <- "lyapunov"

    # Output
    return(LE)
  } else {
    stop("'x' should be a matrix or data frame containing the partial derivatives of the jacobian")
  }
}


