################################################################################
###### DChaos: an R package for detecting chaotic signals in time series #######
################################################################################
# FUNCTION:         CHAOTIC TIME SERIES:
#   logistic.ts       Simulates series from the Logistic equation (discrete time)
#   henon.ts          Simulates series from the Hénon map (discrete time)
#   rossler.ts        Simulates series from the Rössler system (continuous time)
#   lorenz.ts         Simulates series from the Lorenz system (continuous time)
################################################################################
#' Simulation of time series from the Logistic equation
#' @name logistic.ts
#' @aliases logistic.ts
#' @description
#' This function simulates time series from the Logistic equation considering the argument set selected by the user. The initial condition of each time serie is a random number between 0 and 1.
#' @param u.min a non-negative integer denoting a lower bound for the parameter \code{u} (Default 1).
#' @param u.max a non-negative integer denoting an upper bound for the parameter \code{u} (Default 4).
#' @param sample a non-negative integer denoting the length of each time serie (Default 1000).
#' @param transient a non-negative integer denoting the number of observations that will be discarded to ensure that the values of each time serie are in the attractor (Default 100).
#' @param B a non-negative integer denoting the number of series that will be generated for several \code{u}-parameter values. The number of simulated series must be at least 100 (Default 100).
#' @param doplot a logical value denoting if you want to draw a plot \code{TRUE} or not \code{FALSE}. If it is \code{TRUE} shows four graphs: the evolution of the temporal trajectories for an initial period, the evolution of these for the whole period, the attractor and its projections on the Cartesian plane and the bifurcation diagram (Default \code{TRUE}).
#' @return A list containing as many items as series we want to simulate \code{B}. Each of them has the following attributes: the value of the parameter \code{u}, the value of the initial condition \code{xo} and a time serie from the iterated Logistic equation.
#' @references May, R.M. 1976 Simple mathematical models with very complicated dynamics. Nature (261):459-467.
#' @examples
#' ## Simulates 100 time series from the logistic equation for
#' ## u-parameter values in which this system exhibits a chaotic
#' ## behaviour:
#' ts<-logistic.ts(u.min=3.57,u.max=4,B=100,doplot=TRUE)
#' @author Julio E. Sandubete, Lorenzo Escot
#' @seealso \code{\link{henon.ts}}, \code{\link{rossler.ts}}, \code{\link{lorenz.ts}}
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom stats runif
#' @importFrom stats as.ts
#' @export logistic.ts
logistic.ts<-function(u.min=1,u.max=4,sample=1000,transient=100,B=100,doplot=TRUE){

  # Checks
  if ((u.min < 1) | (length(u.min) != 1) | (u.min > u.max)){stop("wrong value for the lower bound of the parameter 'u'")}
  if ((u.max > 4) | (length(u.max) != 1) | (u.max < u.min)){stop("wrong value for the upper bound of the parameter 'u'")}
  if (sample <= 0 | (length(sample) != 1)){stop("wrong value for the length of each time serie")}
  if ((transient <= 0) | (transient > sample) | (length(transient) != 1)){stop("wrong value for the number of observations that will be discarded")}
  if (B<100){stop("the number of simulated series must be at least 100")}

  # Settings
  n<-sample+transient
  u<-seq(from=u.min,to=u.max,length.out=B)
  x<-array(dim=c(n,B))

  # Simulates 'B' time series from the iterated logistic equation
  for(j in 1:B){
    for(i in seq(2,n)){
      x[1,j]<-runif(1,0,1)
      x[i,j]<-u[j]*x[i-1,j]*(1-x[i-1,j])
    }
  }

  # Plots
  if (doplot==TRUE){
    par(mfrow=c(2,2))
    plot(x=seq(1:transient),y=x[1:transient,B],xlab = "Number of observations", ylab = "Values",main=paste("Logistic map | Iteration\n","u =",u.max),type="l",col="darkgray")
    plot(x=seq(1:sample),y=x[(transient+1):n,B],xlab = "Number of observations", ylab = "Values",main=paste("Logistic map | Iteration\n","u =",u.max),type="l",col="darkgray")
    plot(x=x[(transient+1):(n-1),B],y=x[(transient+2):n,B],xlab= "Values (t-1)",ylab="Values (t)",main=paste("Logistic map | Attractor\n","u =",u.max),type="n")
    points(x=x[(transient+1):(n-1),B],y=x[(transient+2):n,B],cex=.34,col="indianred1")
    plot(rep(u[1],sample),x[(transient+1):n,1],xlab="u-parameter values",ylab="Values",main=paste("Logistic map | Bifurcation\n",u.min,"< u <",u.max),xlim=c(u.min,u.max),ylim=c(0,1),type="p",cex=.05,col="indianred1")
    for (k in 2:B){
      points(rep(u[k],sample),x[(transient+1):n,k],cex=.05,col="indianred1")
    }
    par(mfrow=c(1,1))
  }

  # Output
  log.out<-vector("list",B)
  names(log.out)<-paste(rep("Logistic",B),c(1:B))
  for (z in 1:B){
    log.out[[z]]<-list("u"=u[z],"xo"=x[1,z],"time.serie"=stats::as.ts(x[(transient+1):n,z]))
  }
  return(log.out)
}

################################################################################
#' Simulation of time series from the Hénon map
#' @name henon.ts
#' @aliases henon.ts
#' @description
#' This function simulates time series from the Hénon map considering the argument set selected by the user. The initial conditions of each time serie are two random numbers between -0.5 and 0.5. Some values for the parameters and initial conditions may lead to an unstable system that will tend to infinity.
#' @param a.min a non-negative integer denoting a lower bound for the parameter \code{a} (Default 0.4).
#' @param a.max a non-negative integer denoting an upper bound for the parameter \code{a} (Default 1.4).
#' @param b.min a non-negative integer denoting a lower bound for the parameter \code{b} (Default 0.1).
#' @param b.max a non-negative integer denoting an upper bound for the parameter \code{b} (Default 0.3).
#' @param sample a non-negative integer denoting the length of each time serie (Default 1000).
#' @param transient a non-negative integer denoting the number of observations that will be discarded to ensure that the values of each time serie are in the attractor (Default 100).
#' @param B a non-negative integer denoting the number of series that will be generated for different values of the parameters \code{a} and \code{b}. The number of simulated series must be at least 100 (Default 100).
#' @param doplot a logical value denoting if you want to draw a plot \code{TRUE} or not \code{FALSE}. If it is \code{TRUE} shows six graphs: the evolution of the temporal trajectories for the whole period, the attractor and its projections on the Cartesian plane and the bifurcation diagram. All of them consider both the 'x-coordinate' and the 'y-coordinate' (Default 'TRUE').
#' @return A list containing as many items as series we want to simulate \code{B}. Each of them has the following attributes: the value of the parameter \code{a}, the value of the parameter \code{b}, the value of the initial condition \code{xo}, the value of the initial condition \code{yo} and a time serie from the iterated Hénon map with two columns corresponding to 'x-coordinate' and 'y-coordinate'.
#' @references Hénon, M. 1976 A two-dimensional mapping with a strange attractor. Communications in Mathematical Physics 50(1):69-77.
#' @examples
#' ## Simulates 100 time series from the Hénon map for different
#' ## values of the parameters a and b in which this system exhibits
#' ## a chaotic behaviour:
#' ts<-henon.ts(a.min=0.7,a.max=1.4,b.min=0.1,b.max=0.3,B=100,doplot=TRUE)
#' @author Julio E. Sandubete, Lorenzo Escot
#' @seealso \code{\link{logistic.ts}}, \code{\link{rossler.ts}}, \code{\link{lorenz.ts}}
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom stats runif
#' @importFrom stats as.ts
#' @export henon.ts
henon.ts<-function(a.min=0.4,a.max=1.4,b.min=0.1,b.max=0.3,sample=1000,transient=100,B=100,doplot=TRUE){

  # Checks
  if ((a.min < -1.5) | (a.min > 1.5) | (length(a.min) != 1) | (a.min > a.max)){stop("wrong value for the lower bound of the parameter 'a'")}
  if ((a.max < -1.5) | (a.max > 1.5) | (length(a.min) != 1) | (a.max < a.min)){stop("wrong value for the upper bound of the parameter 'a'")}
  if ((b.min < -0.4) | (b.min > 0.4) | (length(b.min) != 1) | (b.min > b.max)){stop("wrong value for the lower bound of the parameter 'b'")}
  if ((b.max < -0.4) | (b.max > 0.4) | (length(b.max) != 1) | (b.max < b.min)){stop("wrong value for the upper bound of the parameter 'b'")}
  if (sample <= 0 | (length(sample) != 1)){stop("wrong value for the length of each time serie")}
  if ((transient <= 0) | (transient > sample) | (length(transient) != 1)){stop("wrong value for the number of observations that will be discarded")}
  if (B<100){stop("the number of simulated series must be at least 100")}

  # Settings
  n<-sample+transient
  a<-seq(from=a.min,to=a.max,length.out=B)
  b<-seq(from=b.min,to=b.max,length.out=B)
  x<-array(dim=c(n,B))
  y<-array(dim=c(n,B))

  # Simulates 'B' time series from the iterated Hénon map
  for(j in 1:B){
    for(i in seq(2,n)){
      x[1,j]<-runif(1,-0.5,0.5)
      y[1,j]<-runif(1,-0.5,0.5)
      x[i,j]<-1-a[j]*x[i-1,j]^2+b[j]*y[i-1,j]
      y[i,j]<-x[i-1,j]
    }
  }

  # To avoid an unstable system
  if (x[n,B] == Inf | x[n,B] == -Inf | y[n,B] == Inf | y[n,B] == -Inf){
    stop("should simulate it again due to some values of the parameters or initial conditions have leaded to an unstable system that tend to infinity.")
  }

  # Plots
  if (doplot==TRUE){
    par(mfrow=c(3,2))
    plot(x=seq(1:sample),y=x[(transient+1):n,B],xlab = "Number of observations", ylab = "Values",main=paste("Henon map (x) | Iteration\n","a =",a.max,",","b =",b.max),type="l",col="darkgray")
    plot(x=seq(1:sample),y=y[(transient+1):n,B],xlab = "Number of observations", ylab = "Values",main=paste("Henon map (y) | Iteration\n","a =",a.max,",","b =",b.max),type="l",col="darkgray")
    plot(x=x[(transient+1):n,B],y=y[(transient+1):n,B],xlab="Henon (x-coordinate)",ylab="Henon (y-coordinate)",main=paste("Henon map | Attractor\n","a =",a.max,",","b =",b.max),type="n")
    points(x=x[(transient+1):n,B],y=y[(transient+1):n,B],cex=.34,col="indianred1")
    plot(x=y[(transient+1):n,B],y=x[(transient+1):n,B],xlab="Henon (y-coordinate)",ylab="Henon (x-coordinate)",main=paste("Henon map | Attractor\n","a =",a.max,",","b =",b.max),type="n")
    points(x=y[(transient+1):n,B],y=x[(transient+1):n,B],cex=.34,col="indianred1")
    plot(rep(a[1],sample),x[(transient+1):n,1],xlab="a-parameter values",ylab="Henon (x-coordinate)",main=paste("Henon map | Bifurcation\n",a.min,"< a <",a.max),xlim=c(a.min,a.max),ylim=c(-1.5,1.5),type="p",cex=.05,col="indianred1")
    for (k in 2:B){
      points(rep(a[k],sample),x[(transient+1):n,k],cex=.05,col="indianred1")
    }
    plot(rep(b[1],sample),y[(transient+1):n,1],xlab="b-parameter values",ylab="Henon (y-coordinate)",main=paste("Henon map | Bifurcation\n",b.min,"< b <",b.max),xlim=c(b.min,b.max),ylim=c(-0.4,0.4),type="p",cex=.05,col="indianred1")
    for (v in 2:B){
      points(rep(b[v],sample),y[(transient+1):n,v],cex=.05,col="indianred1")
    }
    par(mfrow=c(1,1))
  }

  # Output
  hen.out<-vector("list",B)
  names(hen.out)<-paste(rep("Henon",B),c(1:B))
  for (z in 1:B){
    hen.out[[z]]<-list("a"=a[z],"b"=b[z],"xo"=x[1,z],"yo"=y[1,z],"time.serie"=list("Xt"=stats::as.ts(x[(transient+1):n,z]),"Yt"=stats::as.ts(y[(transient+1):n,z])))
  }
  return(hen.out)
}

################################################################################
#' Simulation of time series from the Rössler system
#' @name rossler.ts
#' @aliases rossler.ts
#' @description
#' This function simulates time series from the Rössler system considering the argument set selected by the user. Some values for the parameters and initial conditions may lead to an unstable system that will tend to infinity.
#' @param a.min a non-negative integer denoting a lower bound for the parameter \code{a} (Default 0).
#' @param a.max a non-negative integer denoting an upper bound for the parameter \code{a} (Default 0.2).
#' @param b.min a non-negative integer denoting a lower bound for the parameter \code{b} (Default 0).
#' @param b.max a non-negative integer denoting an upper bound for the parameter \code{b} (Default 0.2).
#' @param c.min a non-negative integer denoting a lower bound for the parameter \code{c} (Default 4).
#' @param c.max a non-negative integer denoting an upper bound for the parameter \code{c} (Default 5.7).
#' @param xo.min a non-negative integer denoting a lower bound for the initial condition \code{xo} (Default -2).
#' @param xo.max a non-negative integer denoting an upper bound for the initial condition \code{xo} (Default 2).
#' @param yo.min a non-negative integer denoting a lower bound for the initial condition \code{yo} (Default -10).
#' @param yo.max a non-negative integer denoting an upper bound for the initial condition \code{yo} (Default 10).
#' @param zo.min a non-negative integer denoting a lower bound for the initial condition \code{zo} (Default -0.2).
#' @param zo.max a non-negative integer denoting an upper bound for the initial condition \code{zo} (Default 0.2).
#' @param time a numeric vector denoting the time-lapse and the time-step (Default 'time-lapse' equal to 10001 with a 'time-step' of 0.01 seconds).
#' @param transient a non-negative integer denoting the number of observations that will be discarded to ensure that the values of each time serie are in the attractor (Default 1000).
#' @param B a non-negative integer denoting the number of series that will be generated for different values of parameters \code{a}, \code{b} and \code{c}. The number of simulated series must be at least 100 (Default 100).
#' @param doplot a logical value denoting if you want to draw a plot \code{TRUE} or not \code{FALSE}. If it is \code{TRUE} shows six graphs: the evolution of the temporal trajectories for the whole period, the attractor and its projections on the Cartesian plane. All of them consider the 'x-coordinate', the 'y-coordinate' and the 'z-coordinate' (Default \code{TRUE}).
#' @return A list containing as many items as series we want to simulate \code{B}. Each of them has the following attributes: the value of the parameter \code{a}, the value of the parameter \code{b}, the value of the parameter \code{c}, the value of the initial condition \code{xo}, the value of the initial condition \code{yo}, the value of the initial condition \code{zo} and a time serie from the iterated Rössler system with three columns corresponding to 'x-coordinate', 'y-coordinate' and 'z-coordinate'.
#' @references Rössler, O. 1976 An equation for continuous chaos. Physics Letters A 57(5):397-398.
#' @examples
#' ## Simulates 100 time series from the Rössler system for different
#' ## values of the parameters a, b and c in which this system exhibits
#' ## a chaotic behaviour:
#' ts<-rossler.ts(a.min=0.2,a.max=0.2,b.min=0.2,b.max=0.2,c.min=5.7,c.max=5.7,
#'                time=seq(0,10,0.01),transient=100,B=100,doplot=TRUE)
#' @author Julio E. Sandubete, Lorenzo Escot
#' @seealso \code{\link{logistic.ts}}, \code{\link{henon.ts}}, \code{\link{lorenz.ts}}
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom stats runif
#' @export rossler.ts
rossler.ts<-function(a.min=0.1,a.max=0.2,b.min=0.1,b.max=0.2,c.min=4,c.max=5.7,xo.min=-2,xo.max=2,yo.min=-10,yo.max=10,zo.min=-0.2,zo.max=0.2,time=seq(0,100,0.01),transient=1000,B=100,doplot=TRUE){

  # Checks
  if ((a.min < 0.1) | (a.min > 10) | (length(a.min) != 1) | (a.min > a.max)){stop("wrong value for the lower bound of the parameter 'a'")}
  if ((a.max < 0.1) | (a.max > 10) | (length(a.max) != 1) | (a.max < a.min)){stop("wrong value for the upper bound of the parameter 'a'")}
  if ((b.min < 0.1) | (b.min > 10) | (length(b.min) != 1) | (b.min > b.max)){stop("wrong value for the lower bound of the parameter 'b'")}
  if ((b.max < 0.1) | (b.max > 10) | (length(b.max) != 1) | (b.max < b.min)){stop("wrong value for the upper bound of the parameter 'b'")}
  if ((c.min < 1) | (c.min > 100) | (length(c.min) != 1) | (c.min > c.max)){stop("wrong value for the lower bound of the parameter 'c'")}
  if ((c.max < 1) | (c.max > 100) | (length(c.max) != 1) | (c.max < c.min)){stop("wrong value for the upper bound of the parameter 'c'")}
  if ((xo.min < -10) | (xo.min > 10) | (length(xo.min) != 1) | (xo.min > xo.max)){stop("wrong value for the lower bound of the initial condition 'xo'")}
  if ((xo.max < -10) | (xo.max > 10) | (length(xo.max) != 1) | (xo.max < xo.min)){stop("wrong value for the upper bound of the initial condition 'xo'")}
  if ((yo.min < -10) | (yo.min > 10) | (length(yo.min) != 1) | (yo.min > yo.max)){stop("wrong value for the lower bound of the initial condition 'yo'")}
  if ((yo.max < -10) | (yo.max > 10) | (length(yo.max) != 1) | (yo.max < yo.min)){stop("wrong value for the upper bound of the initial condition 'yo'")}
  if ((zo.min < -10) | (zo.min > 10) | (length(zo.min) != 1) | (zo.min > zo.max)){stop("wrong value for the lower bound of the initial condition 'zo'")}
  if ((zo.max < -10) | (zo.max > 10) | (length(zo.max) != 1) | (zo.max < zo.min)){stop("wrong value for the upper bound of the initial condition 'zo'")}
  if (!is.vector(time)){stop("'time' should be a numeric vector denoting the time-lapse and the time-step")}
  if ((transient <= 0) | (transient > length(time)) | (length(transient) != 1)){stop("wrong value for the number of observations that will be discarded")}
  if (B<100){stop("the number of simulated series must be at least 100")}

  # Settings
  n<-length(time)
  a<-seq(from=a.min,to=a.max,length.out=B)
  b<-seq(from=b.min,to=b.max,length.out=B)
  c<-seq(from=c.min,to=c.max,length.out=B)
  parameter<-cbind(a,b,c)
  xo<-seq(from=xo.min,to=xo.max,length.out=B)
  yo<-seq(from=yo.min,to=yo.max,length.out=B)
  zo<-seq(from=zo.min,to=zo.max,length.out=B)
  setting<-cbind(xo,yo,zo)

  # Simulates 'B' time series from the iterated Rössler system
  ros.equ<-function(parameter,setting,time){

    # Settings
    a<-parameter[[1]]
    b<-parameter[[2]]
    c<-parameter[[3]]
    xo<-setting[[1]]
    yo<-setting[[2]]
    zo<-setting[[3]]

    # Rössler system
    dx=-yo-zo
    dy=xo+a*yo
    dz=b+zo*(xo-c)

    #Output
    return(c(dx,dy,dz))
  }
  ros.out<-vector("list",B)
  names(ros.out)<-paste(rep("Rossler",B),c(1:B))
  for(j in 1:B){
    ros.sim<-rk(fun=ros.equ,parameter=parameter[j,],setting=setting[j,],time=time)
    ros.out[[j]]<-list("a"=parameter[j,][[1]],"b"=parameter[j,][[2]],"c"=parameter[j,][[3]],
                       "xo"=setting[j,][[1]],"yo"=setting[j,][[2]],"zo"=setting[j,][[3]],
                       "time.serie"=list("Xt"=stats::as.ts(ros.sim$Xt[-c(1:transient)]),"Yt"=stats::as.ts(ros.sim$Yt[-c(1:transient)]),"Zt"=stats::as.ts(ros.sim$Zt[-c(1:transient)])))
  }

  # To avoid an unstable system
  if (is.nan(ros.out[[B]]$time.serie$Xt[n-transient])){stop("should simulate it again due to some values of the parameters or initial conditions have leaded to an unstable system that tend to infinity.")}
  if (is.nan(ros.out[[B]]$time.serie$Yt[n-transient])){stop("should simulate it again due to some values of the parameters or initial conditions have leaded to an unstable system that tend to infinity.")}
  if (is.nan(ros.out[[B]]$time.serie$Zt[n-transient])){stop("should simulate it again due to some values of the parameters or initial conditions have leaded to an unstable system that tend to infinity.")}

  # Plots
  if (doplot==TRUE){
    par(mfrow=c(2,3))
    plot(x=seq(1:(n-transient)),y=ros.out[[B]]$time.serie$Xt,xlab = "Number of observations", ylab = "Values",main=paste("Rossler system (x) | Iteration\n","a =",a.max,",","b =",b.max,",","c =",c.max),type="l",col="darkgray")
    plot(x=seq(1:(n-transient)),y=ros.out[[B]]$time.serie$Yt,xlab = "Number of observations", ylab = "Values",main=paste("Rossler system (y) | Iteration\n","a =",a.max,",","b =",b.max,",","c =",c.max),type="l",col="darkgray")
    plot(x=seq(1:(n-transient)),y=ros.out[[B]]$time.serie$Zt,xlab = "Number of observations", ylab = "Values",main=paste("Rossler system (z) | Iteration\n","a =",a.max,",","b =",b.max,",","c =",c.max),type="l",col="darkgray")
    plot(x=ros.out[[B]]$time.serie$Xt,y=ros.out[[B]]$time.serie$Yt,xlab="Rossler (x-coordinate)",ylab="Rossler (y-coordinate)",main=paste("Rossler system | Attractor\n","a =",a.max,",","b =",b.max,",","c =",c.max),type="n")
    points(x=ros.out[[B]]$time.serie$Xt,y=ros.out[[B]]$time.serie$Yt,cex=.34,col="indianred1")
    plot(x=ros.out[[B]]$time.serie$Yt,y=ros.out[[B]]$time.serie$Zt,xlab="Rossler (y-coordinate)",ylab="Rossler (z-coordinate)",main=paste("Rossler system | Attractor\n","a =",a.max,",","b =",b.max,",","c =",c.max),type="n")
    points(x=ros.out[[B]]$time.serie$Yt,y=ros.out[[B]]$time.serie$Zt,cex=.34,col="indianred1")
    plot(x=ros.out[[B]]$time.serie$Zt,y=ros.out[[B]]$time.serie$Xt,xlab="Rossler (z-coordinate)",ylab="Rossler (x-coordinate)",main=paste("Rossler system | Attractor\n","a =",a.max,",","b =",b.max,",","c =",c.max),type="n")
    points(x=ros.out[[B]]$time.serie$Zt,y=ros.out[[B]]$time.serie$Xt,cex=.34,col="indianred1")
    par(mfrow=c(1,1))
  }

  # Output
  return(ros.out)
}

################################################################################
#' Simulation of time series from the Lorenz system
#' @name lorenz.ts
#' @aliases lorenz.ts
#' @description
#' This function simulates time series from the Lorenz system considering the argument set selected by the user. Some values for the parameters and initial conditions may lead to an unstable system that will tend to infinity.
#' @param sigma.min a non-negative integer denoting a lower bound for the parameter \code{sigma} (Default 8).
#' @param sigma.max a non-negative integer denoting an upper bound for the parameter \code{sigma} (Default 10).
#' @param rho.min a non-negative integer denoting a lower bound for the parameter \code{rho} (Default 25).
#' @param rho.max a non-negative integer denoting an upper bound for the parameter \code{rho} (Default 27).
#' @param beta.min a non-negative integer denoting a lower bound for the parameter \code{beta} (Default 1).
#' @param beta.max a non-negative integer denoting an upper bound for the parameter \code{beta} (Default 2.67).
#' @param xo.min a non-negative integer denoting a lower bound for the initial condition \code{xo} (Default -14).
#' @param xo.max a non-negative integer denoting an upper bound for the initial condition \code{xo} (Default -10).
#' @param yo.min a non-negative integer denoting a lower bound for the initial condition \code{yo} (Default -13).
#' @param yo.max a non-negative integer denoting an upper bound for the initial condition \code{yo} (Default -10).
#' @param zo.min a non-negative integer denoting a lower bound for the initial condition \code{zo} (Default 3).
#' @param zo.max a non-negative integer denoting an upper bound for the initial condition \code{zo} (Default 10).
#' @param time a numeric vector denoting the time-lapse and the time-step (Default 'time-lapse' equal to 10001 with a 'time-step' of 0.01 seconds).
#' @param transient a non-negative integer denoting the number of observations that will be discarded to ensure that the values of each time serie are in the attractor (Default 1000).
#' @param B a non-negative integer denoting the number of series that will be generated for different values of parameters \code{sigma}, \code{rho} and \code{beta}. The number of simulated series must be at least 100 (Default 100).
#' @param doplot a logical value denoting if you want to draw a plot \code{TRUE} or not \code{FALSE}. If it is \code{TRUE} shows six graphs: the evolution of the temporal trajectories for the whole period, the attractor and its projections on the Cartesian plane. All of them consider the 'x-coordinate', the 'y-coordinate' and the 'z-coordinate' (Default \code{TRUE}).
#' @return A list containing as many items as series we want to simulate \code{B}. Each of them has the following attributes: the value of the parameter \code{sigma}, the value of the parameter \code{rho}, the value of the parameter \code{beta}, the value of the initial condition \code{xo}, the value of the initial condition \code{yo}, the value of the initial condition \code{zo} and a time serie from the iterated Lorenz system with three columns corresponding to 'x-coordinate', 'y-coordinate' and 'z-coordinate'.
#' @references Lorenz, E. 1963 Deterministic nonperiodic flow. Journal of the Atmospheric Sciences 20(2):130-141.
#' @examples
#' ## Simulates 100 time series from the Lorenz system for different
#' ## values of the parameters sigma, rho and beta in which this system
#' ## exhibits a chaotic behaviour:
#' ts<-lorenz.ts(sigma.min=10,sigma.max=10,rho.min=27,rho.max=27,beta.min=2.67,
#'               beta.max=2.67,time=seq(0,10,0.01),transient=100,B=100, doplot=TRUE)
#' @author Julio E. Sandubete, Lorenzo Escot
#' @seealso \code{\link{logistic.ts}}, \code{\link{henon.ts}}, \code{\link{rossler.ts}}
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics mtext
#' @importFrom stats runif
#' @export lorenz.ts
lorenz.ts<-function(sigma.min=8,sigma.max=10,rho.min=25,rho.max=27,beta.min=1,beta.max=2.67,xo.min=-14,xo.max=-10,yo.min=-13,yo.max=-10,zo.min=3,zo.max=10,time=seq(0,100,0.01),transient=1000,B=100,doplot=TRUE){

  # Checks
  if ((sigma.min < 0) | (sigma.min > 100) | (length(sigma.min) != 1) | (sigma.min > sigma.max)){stop("wrong value for the lower bound of the parameter 'sigma'")}
  if ((sigma.max < 0) | (sigma.max > 100) | (length(sigma.max) != 1) | (sigma.max < sigma.min)){stop("wrong value for the upper bound of the parameter 'sigma'")}
  if ((rho.min < 0) | (rho.min > 100) | (length(rho.min) != 1) | (rho.min > rho.max)){stop("wrong value for the lower bound of the parameter 'rho'")}
  if ((rho.max < 0) | (rho.max > 100) | (length(rho.max) != 1) | (rho.max < rho.min)){stop("wrong value for the upper bound of the parameter 'rho'")}
  if ((beta.min < 0) | (beta.min > 100) | (length(beta.min) != 1) | (beta.min > beta.max)){stop("wrong value for the lower bound of the parameter 'beta'")}
  if ((beta.max < 0) | (beta.max > 100) | (length(beta.max) != 1) | (beta.max < beta.min)){stop("wrong value for the upper bound of the parameter 'beta'")}
  if ((xo.min < -100) | (xo.min > 100) | (length(xo.min) != 1) | (xo.min > xo.max)){stop("wrong value for the lower bound of the initial condition 'xo'")}
  if ((xo.max < -100) | (xo.max > 100) | (length(xo.min) != 1) | (xo.max < xo.min)){stop("wrong value for the upper bound of the initial condition 'xo'")}
  if ((yo.min < -100) | (yo.min > 100) | (length(yo.min) != 1) | (yo.min > yo.max)){stop("wrong value for the lower bound of the initial condition 'yo'")}
  if ((yo.max < -100) | (yo.max > 100) | (length(yo.max) != 1) | (yo.max < yo.min)){stop("wrong value for the upper bound of the initial condition 'yo'")}
  if ((zo.min < -100) | (zo.min > 100) | (length(zo.min) != 1) | (zo.min > zo.max)){stop("wrong value for the lower bound of the initial condition 'zo'")}
  if ((zo.max < -100) | (zo.max > 100) | (length(zo.max) != 1) | (zo.max < zo.min)){stop("wrong value for the upper bound of the initial condition 'zo'")}
  if (!is.vector(time)){stop("'time' should be a numeric vector denoting the time-lapse and the time-step")}
  if ((transient <= 0) | (transient > length(time)) | (length(transient) != 1)){stop("wrong value for the number of observations that will be discarded")}
  if (B<100){stop("the number of simulated series must be at least 100")}

  # Settings
  n<-length(time)
  sigma<-seq(from=sigma.min,to=sigma.max,length.out=B)
  rho<-seq(from=rho.min,to=rho.max,length.out=B)
  beta<-seq(from=beta.min,to=beta.max,length.out=B)
  parameter<-cbind(sigma,rho,beta)
  xo<-seq(from=xo.min,to=xo.max,length.out=B)
  yo<-seq(from=yo.min,to=yo.max,length.out=B)
  zo<-seq(from=zo.min,to=zo.max,length.out=B)
  setting<-cbind(xo,yo,zo)

  # Simulates 'B' time series from the iterated Lorenz system
  lor.equ<-function(parameter,setting,time){

    # Settings
    sigma<-parameter[[1]]
    rho<-parameter[[2]]
    beta<-parameter[[3]]
    xo<-setting[[1]]
    yo<-setting[[2]]
    zo<-setting[[3]]

    # Lorenz system
    dx=sigma*(yo-xo)
    dy=xo*(rho-zo)-yo
    dz=xo*yo-beta*zo

    #Output
    return(c(dx,dy,dz))
  }
  lor.out<-vector("list",B)
  names(lor.out)<-paste(rep("Lorenz",B),c(1:B))
  for(j in 1:B){
    lor.sim<-rk(fun=lor.equ,parameter=parameter[j,],setting=setting[j,],time=time)
    lor.out[[j]]<-list("sigma"=parameter[j,][[1]],"rho"=parameter[j,][[2]],"beta"=parameter[j,][[3]],
                       "xo"=setting[j,][[1]],"yo"=setting[j,][[2]],"zo"=setting[j,][[3]],
                       "time.serie"=list("Xt"=stats::as.ts(lor.sim$Xt[-c(1:transient)]),"Yt"=stats::as.ts(lor.sim$Yt[-c(1:transient)]),"Zt"=stats::as.ts(lor.sim$Zt[-c(1:transient)])))
  }

  # To avoid an unstable system
  if (is.nan(lor.out[[B]]$time.serie$Xt[n-transient])){stop("should simulate it again due to some values of the parameters or initial conditions have leaded to an unstable system that tend to infinity.")}
  if (is.nan(lor.out[[B]]$time.serie$Yt[n-transient])){stop("should simulate it again due to some values of the parameters or initial conditions have leaded to an unstable system that tend to infinity.")}
  if (is.nan(lor.out[[B]]$time.serie$Zt[n-transient])){stop("should simulate it again due to some values of the parameters or initial conditions have leaded to an unstable system that tend to infinity.")}

  # Plots
  if (doplot==TRUE){
    par(mfrow=c(2,3))
    plot(x=seq(1:(n-transient)),y=lor.out[[B]]$time.serie$Xt,xlab = "Number of observations",ylab = "Values",main="Lorenz system (x) | Iteration",type="l",col="darkgray")
    graphics::mtext(bquote(list(sigma == .(sigma.max), ~ rho == .(rho.max), ~ beta == .(beta.max))),cex=.75)
    plot(x=seq(1:(n-transient)),y=lor.out[[B]]$time.serie$Yt,xlab = "Number of observations",ylab = "Values",main="Lorenz system (y) | Iteration",type="l",col="darkgray")
    graphics::mtext(bquote(list(sigma == .(sigma.max), ~ rho == .(rho.max), ~ beta == .(beta.max))),cex=.75)
    plot(x=seq(1:(n-transient)),y=lor.out[[B]]$time.serie$Zt,xlab = "Number of observations",ylab = "Values",main="Lorenz system (z) | Iteration",type="l",col="darkgray")
    graphics::mtext(bquote(list(sigma == .(sigma.max), ~ rho == .(rho.max), ~ beta == .(beta.max))),cex=.75)
    plot(x=lor.out[[B]]$time.serie$Xt,y=lor.out[[B]]$time.serie$Yt,xlab="Lorenz (x-coordinate)",ylab="Lorenz (y-coordinate)",main="Lorenz system | Attractor",type="n")
    graphics::mtext(bquote(list(sigma == .(sigma.max), ~ rho == .(rho.max), ~ beta == .(beta.max))),cex=.75)
    points(x=lor.out[[B]]$time.serie$Xt,y=lor.out[[B]]$time.serie$Yt,cex=.34,col="indianred1")
    plot(x=lor.out[[B]]$time.serie$Yt,y=lor.out[[B]]$time.serie$Zt,xlab="Lorenz (y-coordinate)",ylab="Lorenz (z-coordinate)",main="Lorenz system | Attractor",type="n")
    graphics::mtext(bquote(list(sigma == .(sigma.max), ~ rho == .(rho.max), ~ beta == .(beta.max))),cex=.75)
    points(x=lor.out[[B]]$time.serie$Yt,y=lor.out[[B]]$time.serie$Zt,cex=.34,col="indianred1")
    plot(x=lor.out[[B]]$time.serie$Zt,y=lor.out[[B]]$time.serie$Xt,xlab="Lorenz (z-coordinate)",ylab="Lorenz (x-coordinate)",main="Lorenz system | Attractor",type="n")
    graphics::mtext(bquote(list(sigma == .(sigma.max), ~ rho == .(rho.max), ~ beta == .(beta.max))),cex=.75)
    points(x=lor.out[[B]]$time.serie$Zt,y=lor.out[[B]]$time.serie$Xt,cex=.34,col="indianred1")
    par(mfrow=c(1,1))
  }

  # Output
  return(lor.out)
}

################################################################################
# Internal function: approximate solutions of ordinary differential equations.
# We have implemented the classical Runge–Kutta method (RK4). It is used to
# generate time series from both Rössler system and Lorenz system.
rk<-function(fun,parameter,setting,time){

  # Checks
  if (!is.function(fun)){stop("'fun' should be a function containing the Rossler system or the Lorenz system")}
  if (!is.vector(parameter)){stop("'parameter' should be a numeric vector denoting the parameter set")}
  if (!is.vector(setting)){stop("'setting' should be a numeric vector denoting the setting set")}
  if (!is.vector(time)){stop("'time' should be a numeric vector denoting the time-lapse and the time-step")}

  # Settings
  n<-length(time)
  out<-matrix(ncol=length(setting), nrow=n)
  out[1,]<-setting
  h<-time[2]-time[1]

  # Runge-Kutta method (RK4)
  for (i in 2:n){
    k1<-h*fun(parameter,out[i-1,],time[i-1])
    k2<-h*fun(parameter,out[i-1,]+k1/2,time[i-1]+h/2)
    k3<-h*fun(parameter,out[i-1,]+k2/2,time[[i-1]]+h/2)
    k4<-h*fun(parameter,out[i-1,]+k3,time[[i-1]]+h)
    out[i,]<-out[i-1,]+(1/6)*(k1+2*k2+2*k3+k4)
  }

  # Output
  return(list("Xt"=stats::as.ts(out[,1]),"Yt"=stats::as.ts(out[,2]),"Zt"=stats::as.ts(out[,3])))
}










