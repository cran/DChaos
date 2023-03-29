################################################################################
############ DChaos: An R Package for Chaotic Time Series Analysis #############
################################################################################
# FUNCTION:         DESCRIPTION:
#   lyapunov.max      Estimates the largest Lyapunov Exponent (NORMA-2)
################################################################################
#' Estimates the largest Lyapunov exponent
#' @name lyapunov.max
#' @aliases lyapunov.max
#' @description
#' This function estimates the largest Lyapunov exponent through the Norma-2 procedure based on the partial derivatives computed by the \code{jacobian.net} function.
#' @param data should be a \code{jacobian} object containing the partial derivatives computed by the \code{jacobian.net} function.
#' @param blocking a character denoting the blocking method chosen for figuring out the largest Lyapunov exponent through the Norma-2 procedure. Available options are \code{FULL} if the user considers the full sample, \code{NOVER} if the user considers the non-overlapping sample, \code{EQS} if the user considers the equally spaced sample or \code{BOOT} if the user considers the bootstrap sample (Default \code{BOOT}).
#' @param B a non-negative integer denoting the number of bootstrap iterations (Default 1000).
#' @param doplot a logical value denoting if the user wants to draw a plot \code{TRUE} or not \code{FALSE}. If it is \code{TRUE} the evolution of the Lyapunov exponent values are represented for the whole period considering the blocking method chosen by the user. It shows as many graphs as embedding dimensions have been considered (Default \code{TRUE}).
#' @return This function returns several objects considering the parameter set selected by the user. The largest Lyapunov exponent considering the Norma-2 procedure by each blocking method are estimated. It also contains some useful information about the estimated jacobian, the best-fitted feed-forward single hidden layer neural net model, the best set of weights found, the fitted values, the residuals obtained, the best embedding parameters set chosen, the sample size or the block length considered by each blocking method. This function provides the standard error, the z test value and the p-value for testing the null hypothesis \eqn{H0: \lambda_k > 0 for k = 1} (largest). Reject the null hypothesis ${H_0}$ means lack of chaotic behaviour. That is, the data-generating process does not have a chaotic attractor because of it does not show the property of sensitivity to initial conditions.
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
#' ## Provides the largest Lyapunov exponent by the Norma-2 procedure considering the
#' ## bootstrap blocking method from the best-fitted neural net model and the partial
#' ## derivatives showed in the jacobian.net example.
#' ## jacobian <- DChaos::jacobian.net(data=ts, m=3:3, lag=1:1, timelapse="FIXED", h=2:10)
#' ## summary(jacobian)
#' ## exponent <- DChaos::lyapunov.max(data=jacobian, blocking="BOOT", B=100, doplot=FALSE)
#' ## summary(exponent)
#' @author Julio E. Sandubete, Lorenzo Escot
#' @importFrom pracma normest
#' @importFrom pracma zeros
#' @importFrom sandwich lrvar
#' @importFrom stats pnorm
#' @importFrom graphics abline
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom stats median
#' @importFrom stats IQR
#' @export lyapunov.max
lyapunov.max     <- function(data, blocking=c("BOOT","NOVER","EQS","FULL"), B=1000, doplot=TRUE){

  # Checks
  if (is.null(data$jacobian)){stop("'data' should be a jacobian object")}
  if (is.matrix(data$jacobian) | is.data.frame(data$jacobian)){
  if (is.null(blocking)){stop("'blocking' should be 'BOOT', 'NOVER', 'EQS' or 'FULL'")}
  if (B<1){stop("wrong value of bootstrap iterations")}

    # Settings
    x <- data$jacobian
    N <- nrow(x)
    m <- ncol(x)
    blocking = match.arg(blocking)

    # Estimates the Lyapunov Exponent: Bootstrap sample
    if (blocking == "BOOT"){
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
        lvmax_mean<-mean(lvmax[M,abs(IQR(lvmax[M,]))<1.5])
        sdlvmax_mean<-mean(varmax[1,abs(IQR(varmax[1,]))<1.5])
        Ztestmax_mean<-lvmax_mean/(sdlvmax_mean/sqrt(M))
        p.value.max_mean<-pnorm(Ztestmax_mean)

        # Median value
        lvmax_median<-median(lvmax[M,])
        sdlvmax_median<-median(varmax[1,])
        Ztestmax_median<-lvmax_median/(sdlvmax_median/sqrt(M))
        p.value.max_median<-pnorm(Ztestmax_median)

        # Output
        LE<-list("estimator"=c("Largest Lyapunov exponent"),"procedure"=c("Norma-2 by bootstrap blocking method"),
                 "exponent.mean"=matrix(c(lvmax_mean,sdlvmax_mean,Ztestmax_mean,p.value.max_mean),nrow=1,ncol=4, byrow=F, dimnames=list("Exponent",c("Estimate","Std. Error","z value","Pr(>|z|)"))),
                 "exponent.median"=matrix(c(lvmax_median,sdlvmax_median,Ztestmax_median,p.value.max_median),nrow=1,ncol=4, byrow=F, dimnames=list("Exponent",c("Estimate","Std. Error","z value","Pr(>|z|)"))),
                 "sample"=N,"block.length"=M,"no.block"=B
        )
        return(LE)
      }
      LE = lyap_max_boot(x,B=B,doplot=doplot)
    }

    # Estimates the Lyapunov Exponent: Non-overlapping sample
    if (blocking == "NOVER"){
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
        LE<-list("estimator"=c("Largest Lyapunov exponent"),"procedure"=c("Norma-2 by non-overlapping blocking method"),
                 "exponent.mean"=matrix(c(lvmax_mean,sdlvmax_mean,Ztestmax_mean,p.value.max_mean),nrow=1,ncol=4, byrow=F, dimnames=list("Exponent",c("Estimate","Std. Error","z value","Pr(>|z|)"))),
                 "exponent.median"=matrix(c(lvmax_median,sdlvmax_median,Ztestmax_median,p.value.max_median),nrow=1,ncol=4, byrow=F, dimnames=list("Exponent",c("Estimate","Std. Error","z value","Pr(>|z|)"))),
                 "sample"=N,"block.length"=M,"no.block"=B
        )
        return(LE)
      }
      LE = lyap_max_over(x,doplot=doplot)
    }

    # Estimates the Lyapunov Exponent: Equally spaced sample
    if (blocking == "EQS"){
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
        LE<-list("estimator"=c("Largest Lyapunov exponent"),"procedure"=c("Norma-2 by equally spaced blocking method"),
                 "exponent.mean"=matrix(c(lvmax_mean,sdlvmax_mean,Ztestmax_mean,p.value.max_mean),nrow=1,ncol=4, byrow=F, dimnames=list("Exponent",c("Estimate","Std. Error","z value","Pr(>|z|)"))),
                 "exponent.median"=matrix(c(lvmax_median,sdlvmax_median,Ztestmax_median,p.value.max_median),nrow=1,ncol=4, byrow=F, dimnames=list("Exponent",c("Estimate","Std. Error","z value","Pr(>|z|)"))),
                 "sample"=N,"block.length"=M,"no.block"=B
        )
        return(LE)
      }
      LE = lyap_max_ess(x,doplot=doplot)
    }

    # Estimates the Lyapunov Exponent: Full sample
    if (blocking == "FULL"){
      lyap_max_full<-function(x,doplot=TRUE){

        # Lyapunov exponents
        lvmax<-matrix(0,N,1)
        J<-matrix(0,m,m)
        TM<-diag(1,m)
        v<-matrix(0,m,1)
        v[1,1]<-1
        for (i in 1:N){
          J<-rbind(as.numeric(x[i,1:m]),diag(1,m-1,m))
          TM0<-TM
          TM<-J%*%TM
          norsal<-try(pracma::normest(TM%*%v),silent = FALSE)
          if('try-error' %in% class(norsal)){
            TM<-TM0
            lvmax[i]<-lvmax[i-1]
            warning("Lyapunov exponent tends to infinity. Not convergence \n")
            next
          } else {
            lvmax[i]<-(1/i)*log(norsal)
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
        LE<-list("estimator"=c("Largest Lyapunov exponent"),"procedure"=c("Norma-2 by full sample method"),
                 "exponent"=matrix(c(lvmax[N],varmax,Ztestmax,p.value.max), nrow=1,ncol=4, byrow=F, dimnames=list("Exponent",c("Estimate","Std. Error","z value","Pr(>|z|)"))),
                 "sample"=N,"block.length"=N,"no.block"=1
        )
        return(LE)
      }
      LE = lyap_max_full(x,doplot=doplot)
    }

    # Class definition
    LE <- c(data,LE,nprint=0)
    class(LE) <- "lyapunov"

    # Output
    return(LE)
  } else {
    stop("'data' should be a jacobian object")
  }
}
