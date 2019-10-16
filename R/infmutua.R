################################################################################
###### DChaos: an R package for detecting chaotic signals in time series #######
################################################################################
# FUNCTION:         DESCRIPTION:
#   infmutua         Estimate the Average Mutual Information function
################################################################################
#' Estimation of the Average Mutual Information function
#' @name infmutua
#' @aliases infmutua
#' @description
#' This function estimates the Average Mutual Information function considering the argument set selected by the user.
#' @param x a numeric vector or time serie.
#' @param partitions a non-negative integer denoting the number of grouping of the set's elements into non-empty subsets, in such a way that every element is included in exactly one subset.
#' @param lag.max a non-negative integer denoting an upper bound for the reconstruction delay (Default 20).
#' @param doplot a logical value denoting if you want to draw a plot \code{TRUE} or not \code{FALSE}.
#' @return The optimum lag which corresponds with the first minimum of the Average Mutual Information function.
#' @references Fraser, A.M., Swinney, H.L. 1986 Independent coordinates for strange attractors from mutual information. Physical Review A 33(2):1134.
#' @examples
#' ## The first minimum of the average mutual information
#' ## function is showed by simulating a time series from
#' ## the logistic equation.
#' data<-logistic.ts(u.min=4,u.max=4,B=100,doplot=FALSE)
#' ts<-data$`Logistic 100`$time.serie
#' lag.opt<-infmutua(ts,lag.max=10)
#' show(lag.opt$MutualInf)
#' show(lag.opt$FirstMin)
#' @author Julio E. Sandubete, Lorenzo Escot
#' @importFrom entropy discretize2d
#' @importFrom entropy mi.plugin
#' @importFrom entropy chi2indep.plugin
#' @importFrom stats chisq.test
#' @export infmutua
infmutua<-function(x,partitions=ceiling(1.5+log(length(x))/log(2)),lag.max=20,doplot=TRUE){

  # Settings
  datos<-matrix(0,ncol=4,nrow=lag.max+1, dimnames = list(NULL,c("lag", "InfMutua", "InfMutChisq/2", "Chisq.test_p-value")))

  # Average Mutual Information
  for (lag in 0:lag.max){
    datos[lag+1,1]<-lag
    xc<-entropy::discretize2d(x[(lag+1):length(x)],x[1:(length(x)-lag)],numBins1 = partitions, numBins2 =partitions)
    datos[lag+1,2]<-entropy::mi.plugin(xc)
    datos[lag+1,3]<-entropy::chi2indep.plugin(xc)*0.5
    datos[lag+1,4]<-round(chisq.test(xc)$p.value,6)
  }

  # Plot
  if(doplot==TRUE)
    {plot(datos[,2]~datos[,1], type="b", xlab="lags", ylab="Average Mutual Information values")}

  # First min
  minindex=1
  for (i in 1:(nrow(datos)-1)){
    if (datos[i,2]<datos[i+1,2]){
      minindex= i
      break
    }
  }
  if (minindex==1) {
    firstminlag=1
  } else {
    firstminlag=minindex-1
  }

  cat(paste0("first minimum: ",firstminlag))

  return(list(MutualInf=datos, FirstMin=firstminlag))
}

