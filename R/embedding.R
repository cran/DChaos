################################################################################
###### DChaos: an R package for detecting chaotic signals in time series #######
################################################################################
# FUNCTION:         DESCRIPTION:
#   embedding         Generates the embedding vectors using the method of delays
################################################################################
#' Construction of embedding vectors using the method of delays
#' @name embedding
#' @aliases embedding
#' @description
#' This function generates both the uniform and non-uniform embedding vectors set from an univariate time serie considering the argument set selected by the user.
#' @param x a numeric vector, time serie, data frame or matrix depending on the method selected in \code{timelapse}.
#' @param m a non-negative integer denoting the embedding dimension (Default 3).
#' @param lag a non-negative integer denoting the reconstruction delay (Default 1).
#' @param timelapse a character denoting if you consider that the observations are sampled at uniform time intervals \code{FIXED} or with a variable time-lapse between each observation \code{VARIABLE} (Default \code{FIXED}).
#' @details If \code{FIXED} has been selected \code{x} must be a numeric vector or time serie. Otherwise \code{VARIABLE} has to be specified. In this case \code{x} must be a data frame or matrix with two columns. First, the date with the following format \code{YMD H:M:OS3} considering milliseconds e.g., 20190407 00:00:03.347. If you don't consider milliseconds you must put .000 after the seconds. It should be an object of class \code{Factor}. Second, the univariate time serie as a sequence of numerical values.
#' @return A data frame with the uniform or non-uniform embedding vectors by columns from an univariate time serie considering the parameter set selected by the user.
#' @references Ruelle, D., Takens, F. 1971 On the nature of turbulence. Communications in Mathematical Physics 20(3):167-192.
#' @references Packard, N.H., Crutchfield, J.P., Farmer, J.D., Shaw, R.S. 1980 Geometry from a time serie. Physical Review Letters 45:712-716.
#' @references Takens, F. 1981 Detecting strange attractors in turbulence. Springer Berlin Heidelberg.
#' @references Sauer, T., Yorke, J.A., Casdagli, M. 1991 Embedology. Journal of Statistical Physics 65(3):579-616.
#' @references Huke, J.P., Broomhead, D.S. 2007 Embedding theorems for non-uniformly sampled dynamical systems. Nonlinearity 20(9):205-244.
#' @examples
#' ## The first ten values corresponding to the uniform embedding
#' ## vectors set for m=5 and lag=1 are showed by simulating
#' ## time series from the logistic equation.
#' data<-logistic.ts(u.min=4,u.max=4,B=100,doplot=FALSE)
#' ts<-data$`Logistic 100`$time.serie
#' embed<-embedding(ts,m=5,lag=1,timelapse="FIXED")
#' show(head(embed, 10))
#' @author Julio E. Sandubete, Lorenzo Escot
#' @seealso \code{\link{jacobi}}
#' @importFrom stats na.omit
#' @importFrom xts cbind.xts
#' @importFrom xts is.xts
#' @importFrom outliers scores
#' @importFrom zoo index
#' @importFrom zoo coredata
#' @export embedding
embedding<-function(x,m=3,lag=1,timelapse=c("FIXED","VARIABLE")){

  # Checks
  if (is.null(x)){stop("'x' should be a numeric vector, time serie, data frame or matrix depending on the method selected in 'timelapse'")}
  if (m < 1){stop("wrong value for the embedding dimension")}
  if (lag < 1){stop("wrong value for the reconstruction delay")}
  if (is.null(timelapse)){stop("'timelapse' should be 'FIXED' or 'VARIABLE'")}
  timelapse = match.arg(timelapse)

  # Method 1: Uniform embedding vectors
  embed.cte<-function(x,m=3,lag=1){

    if (xts::is.xts(x)){stop("'x' should be a numeric vector, time serie, data frame or matrix")}
    if (is.data.frame(x) | is.matrix(x)){
      xout<-as.vector(x[,1])
    } else {
      xout<-x
    }
    if (m == 1){
      xout<-xout[-(1:lag)]
      return(xout)
      stop()
    } else {
      xout<-xout[-(1:((m-1)*lag))]
      l<-length(xout)-1
      for (i in 1:(m-1)){
        xout<-cbind(xout,x[(((m-i)*lag)-(lag-1)):((m-i)*lag+l-(lag-1))])
      }
      xout<-as.data.frame(xout)
      colnames(xout)[1]<-"Xt"
      name.c<-paste(make.unique(rep("Xt-", m), sep = ""),"lag",sep="")
      name.r<-make.unique(rep("t=", (nrow(xout)+1)), sep = "")
      colnames(xout)[c(2:m)]<-name.c[c(2:m)]
      rownames(xout)<-name.r[c(2:(nrow(xout)+1))]
      return(xout)
    }
  }

  # Settings
  if (timelapse == "FIXED")
    xout = embed.cte(x,m=m,lag=lag)

  # Method 2: Non-uniform embedding vectors
  embed.vble<-function(x,m=3,lag=1){

    if (is.data.frame(x) | is.matrix(x)){
      base::options(digits.secs=3)
      x[,1]<-base::as.POSIXct(base::paste0(base::substr(x[,1],1,11),base::substr(x[,1],13,14),base::substr(x[,1],16,17),".",substr(x[,1],19,21)),format="%Y%m%d %H%M%OS")
      x<-xts::xts(x[,2],order.by=x[,1])
      date<-zoo::index(x)
      freq<-diff(date, units = "seconds")
      freq<-c(0,freq)
      x<-xts::cbind.xts(x[,1],freq)
      x<-x[abs(outliers::scores(x[,2], type="z"))<=5]
      x<-xts::cbind.xts(x,cumsum(x[,2]))
      n<-nrow(x)
      tt<-c()
      l<-c()
      z<-c()
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
      xout<-na.omit(xout)
      xout<-as.data.frame(zoo::coredata(xout))
      colnames(xout)[1]<-"Xt"
      name.c<-paste(make.unique(rep("Xt-", m), sep = ""),"lag",sep="")
      name.r<-make.unique(rep("t=", (nrow(xout)+1)), sep = "")
      colnames(xout)[c(2:m)]<-name.c[c(2:m)]
      rownames(xout)<-name.r[c(2:(nrow(xout)+1))]
      return(xout)
    } else {
      stop("'x' should be a data frame or matrix with two columns, see R Documentation")
    }
  }

  # Settings
  if (timelapse == "VARIABLE")
    xout = embed.vble(x,m=m,lag=lag)

  # Output
  return(xout)
}
