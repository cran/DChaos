################################################################################
############ DChaos: An R Package for Chaotic Time Series Analysis #############
################################################################################
# FUNCTION:         DESCRIPTION:
#   embedding         Provides the delayed-coordinate embedding vectors backwards
################################################################################
#' Provides the delayed-coordinate embedding vectors backwards
#' @name embedding
#' @aliases embedding
#' @description
#' This function generates both the uniform and non-uniform embedding vectors backwards using the method of delays from univariate time-series data.
#' @param x a \code{vector}, a time-series object \code{ts} or \code{xts}, a \code{data.frame}, a \code{data.table} or a \code{matrix} depending on the method selected in \code{timelapse}.
#' @param m a non-negative integer denoting the embedding dimension (Default 2).
#' @param lag a non-negative integer denoting the reconstruction delay (Default 1).
#' @param timelapse a character denoting if the time-series data are sampled at uniform time-frequency e.g., 1-month, 1-day, 1-hour, 30-min, 5-min, 1-min and so on \code{FIXED} or non-uniform time-frequency which are not equally spaced in time \code{VARIABLE} (Default \code{FIXED}).
#' @return The uniform or non-uniform delayed-coordinate embedding vectors backwards by columns from an univariate time-series data considering the parameter set selected by the user. If \code{FIXED} has been selected \code{data} must be a \code{vector} or a time-series object \code{ts} or \code{xts}. Otherwise \code{VARIABLE} has to be specified. In this case \code{data} must be a \code{data.frame}, a \code{data.table} or a \code{matrix} with two columns, the date and the univariate time series as a sequence of numerical values, in that order. The date can have the following three classes: \code{POSIXt}, \code{Date} or \code{Factor}. In the latter case the date should come in the following format \code{YMD H:M:OS3} considering milliseconds e.g., 20190407 00:00:03.347. If you don't consider milliseconds you must put .000 after the seconds.
#' @note Note that a key point to create a suitable reconstruction of the state-space is to fix a criteria in order to estimate the embedding parameters. Researchers usually estimate them using heuristic approaches based on prescriptions proposed by e.g., H.D. Abarbanel (1996) or H. Kantz and T. Schreiber (2004). The main drawbacks of these heuristic approaches are the following: they are not intrinsically statistical; their results are not robust; they lead to estimators whose properties are unknown or largely unexplored; they do not take into account the results of any model fit. The alternative proposed by the statistical approach solves those disadvantages. The statistical approach to state-space reconstruction can be viewed as a best subset selection problem within the nonparametric regression context as argued K.-S. Chan and H. Tong (2001). The DChaos package allows the R users to choose between both methods. By default it uses the statistical approach based on model selection procedures instead of heuristic techniques, see \code{netfit} function.
#' @references Ruelle, D., Takens, F. 1971 On the nature of turbulence. Communications in Mathematical Physics 20(3):167-192.
#' @references Takens, F. 1981 Detecting strange attractors in turbulence. Springer Berlin Heidelberg.
#' @references Abarbanel, H.D. 1996 Analysis of observed chaotic data. Springer.
#' @references Cha, K.-S., Tong, H. 2001 Chaos: a statistical perspective. Springer-Verlag.
#' @references Kantz, H., Schreiber, T. 2004 Nonlinear time series analysis, volume 7. Cambridge university press.
#' @references Huke, J.P., Broomhead, D.S. 2007 Embedding theorems for non-uniformly sampled dynamical systems. Nonlinearity 20(9):205-244.
#' @examples
#' ## set.seed(34)
#' ## Simulates time-series data from the Logistic map with chaos
#' ## ts        <- DChaos::logistic.sim(n=1000, a=4)
#' ## show(head(ts, 5))
#'
#' ## Provides the uniform delayed-coordinate embedding vectors (Backward)
#' ## data      <- DChaos::embedding(ts, m=5, lag=2, timelapse="FIXED")
#' ## show(head(data, 5))
#'
#' ## Simulates tick-by-tick data (bid price) for Starbucks company
#' ## ts        <- highfrequency::sbux
#' ## show(head(ts, 5))
#'
#' ## Provides the non-uniform delayed-coordinate embedding vectors (Backward)
#' ## data      <- DChaos::embedding(ts, m=3, lag=4, timelapse="VARIABLE")
#' ## show(head(data, 5))
#' @author Julio E. Sandubete, Lorenzo Escot
#' @importFrom stats na.omit
#' @importFrom xts cbind.xts
#' @importFrom xts is.xts
#' @importFrom zoo index
#' @importFrom zoo coredata
#' @export embedding
embedding<-function(x, m=2, lag=1, timelapse=c("FIXED","VARIABLE")){

  # Checks
  if (is.null(x)){stop("'x' should be a vector, a time-series object ts or xts, a data.frame, a data.table or a matrix depending on the method selected in 'timelapse'")}
  if (m < 1){stop("wrong value for the embedding dimension")}
  if (lag < 1){stop("wrong value for the reconstruction delay")}
  if (is.null(timelapse)){stop("'timelapse' should be 'FIXED' or 'VARIABLE'")}
  timelapse = match.arg(timelapse)

  # Method 1: Uniform delayed-coordinate embedding vectors
  embed.cte<-function(x,m=2,lag=1){

    if (xts::is.xts(x)){stop("'x' should be a vector, a time-series object ts, a data.frame, a data.table or a matrix")}
    if (is.data.frame(x) | is.matrix(x)){
      xout<-as.vector(x[,1])
    } else {
      xout<-x
    }
    if (m == 1){
      xout<-as.data.frame(xout[-(1:lag)])
      colnames(xout)[1]<-"y"
      return(xout)
      stop()
    } else {
      xout<-xout[-(1:((m-1)*lag))]
      l<-length(xout)-1
      for (i in 1:(m-1)){
        xout<-cbind(xout,x[(((m-i)*lag)-(lag-1)):((m-i)*lag+l-(lag-1))])
      }
      xout<-as.data.frame(xout)
      colnames(xout)<-c("y",paste0("x", 1:(m-1)))
      return(xout)
    }
  }

  # Settings
  if (timelapse == "FIXED")
    xout = embed.cte(x,m=m,lag=lag)

  # Method 2: Non-uniform delayed-coordinate embedding vectors
  embed.vble<-function(x,m=2,lag=1){

    if (is.data.frame(x) | is.matrix(x) | is.xts(x)){
      base::options(digits.secs=3)
      if (is.list(x)){
        x<-xts::xts(x[[2]],order.by=x[[1]])
      } else {
        if (is.xts(x)){
          x<-x
        } else {
          x[,1]<-base::as.POSIXct(base::paste0(base::substr(x[,1],1,11),base::substr(x[,1],13,14),base::substr(x[,1],16,17),".",substr(x[,1],19,21)),format="%Y%m%d %H%M%OS")
          x<-xts::xts(x[,2],order.by=x[,1])
        }
      }
      date<-zoo::index(x)
      freq<-diff(date, units = "seconds")
      freq<-c(0,freq)
      x<-xts::cbind.xts(x[,1],freq)
      x<-x[abs(scale(x[,2])) <= 5]
      x<-xts::cbind.xts(x,cumsum(x[,2]))
      colnames(x)<-c("x","difftime","cum.difftime")
      n<-nrow(x)
      tt<-c()
      l<-c()
      z<-c()
      if (m == 1){
        ii<-1
        for(i in 1:n){
          tt[i]<-x[i,3]-(1*lag)
          if (tt[i]<0) l[i]<-0 else {
            while (tt[i]>x[ii,3]) ii<-ii+1
            if (tt[i]==x[ii,3]) l[i]<-x[ii,1] else {
              l[i]<-x[ii-1,1]
            }
          }
        }
        z<-which(tt >= 0)[1]
        xout<-as.data.frame(l[-(1:(z-1))])
        colnames(xout)[1]<-"y"
        return(xout)
        stop()
      } else {
        xt<-matrix(nrow=n,ncol=m-1)
        for(j in 1:(m-1)){
          ii<-c()
          ii<-1
          for(i in 1:n){
            tt[i]<-x[i,3]-(j*lag)
            if (tt[i]<0) l[i]<-0 else {
              while (tt[i]>x[ii,3]) ii<-ii+1
              if (tt[i]==x[ii,3]) l[i]<-x[ii,1] else {
                l[i]<-x[ii-1,1]
              }
            }
          }
          z[j]<-which(tt >= 0)[1]
          xt[,j]<-l
        }
        xout<-cbind(x[-(1:(z[j]-1)),1],xt[-(1:(z[j]-1)),])
        colnames(xout)<-c("y",paste0("x", 1:(m-1)))
        return(xout)
      }
    } else {
      stop("'x' should be a data.table, data.frame, matrix or xts object with two columns, see R Documentation")
    }
  }

  # Settings
  if (timelapse == "VARIABLE")
    xout = embed.vble(x,m=m,lag=lag)

  # Output
  return(xout)
}
