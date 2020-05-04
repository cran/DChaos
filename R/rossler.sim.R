################################################################################
############ DChaos: An R Package for Chaotic Time Series Analysis #############
################################################################################
# FUNCTION:         DESCRIPTION:
#   rossler.sim       Simulates time series from the Rossler system
################################################################################
#' Simulates time-series data from the Rossler system
#' @name rossler.sim
#' @aliases rossler.sim
#' @description
#' This function simulates time-series data from the Rossler system considering the parameter set selected by the user. The initial condition is a random number from the normal distribution with mean equal to 0 and variance equal to 1. Some initial conditions may lead to an unstable system that will tend to infinity.
#' @param a a non-negative integer denoting the value of parameter \code{a} (Default 0.2).
#' @param b a non-negative integer denoting the value of parameter \code{b} (Default 0.2).
#' @param c a non-negative integer denoting the value of parameter \code{c} (Default 5.7).
#' @param s a non-negative integer denoting the variance value of the error term. If \eqn{s=0} gives the standard deterministic map (Default 0).
#' @param x0 a non-negative integer denoting the initial condition of x-coordinate (Default random number from the normal distribution).
#' @param y0 a non-negative integer denoting the initial condition of y-coordinate (Default random number from the normal distribution).
#' @param z0 a non-negative integer denoting the initial condition of z-coordinate (Default random number from the normal distribution).
#' @param time a numeric vector denoting the time-lapse and the time-step (Default time-lapse equal to 10000 with a time-step of 0.01 seconds)
#' @param n.start a non-negative integer denoting the number of observations that will be discarded to ensure that the values are in the attractor (Default 50).
#' @return A time-series data object generated from the Rossler system with or without an additive measurement noise term. This dataset could be useful for researchers interested in the field of chaotic dynamic systems and non-linear time series analysis and professors (and students) who teach (learn) courses related to those topics.
#' @note This function provides also noisy time-series data from the deterministic rossler system adding an additive measurement noise term if \eqn{s>0}. We have added to each time-series data a normal multinomial error term denoted by \eqn{{\varepsilon _t} \sim N\left( {0,s} \right)} with different variance values (\eqn{s}). In this sense we have considered it appropriate to add a measurement noise term because most real-world observed time-series data are usually noise-contaminated signals, characterised by an erratic and persistent volatility in certain periods and there is almost always a source of noise linked to measurement errors in real-world datasets. It has been implemented the classical Runge–Kutta method (RK4) in order to generate time-series data from continuous-time dynamical system as the Rossler system.
#' @references Rössler, O. 1976 An equation for continuous chaos. Physics Letters A 57(5):397-398.
#' @examples
#' ## set.seed(34)
#' ## Simulates time-series data from the deterministic rossler system
#' ## with a chaotic behaviour.
#' ts <- rossler.sim(a=0.2, b=0.2, c=5.7, s=0, time=seq(0,100,0.1))
#' ##
#' ## Simulates time-series data from the deterministic rossler system
#' ## with a non-chaotic behaviour.
#' ts <- rossler.sim(a=0.1, b=0.1, c=7, s=0, time=seq(0,100,0.1))
#' @author Julio E. Sandubete, Lorenzo Escot
#' @importFrom stats rnorm
#' @importFrom stats as.ts
#' @importFrom stats complete.cases
#' @export rossler.sim
rossler.sim<-function(a=0.2, b=0.2, c=5.7, s=0, x0=rnorm(1), y0=rnorm(1), z0=rnorm(1), time=seq(0,100,0.01), n.start = 50){

  # Checks
  if (!is.vector(time)){stop("'time' should be a numeric vector denoting the time-lapse and the time-step")}

  # Settings
  n<-length(time)
  parameter<-c(a,b,c)
  setting<-c(x0,y0,z0)

  # Simulates time series from the Rössler system
  ros.equ<-function(parameter,setting,time){

    # Settings
    a<-parameter[1]
    b<-parameter[2]
    c<-parameter[3]
    xo<-setting[1]
    yo<-setting[2]
    zo<-setting[3]

    # Rössler system
    dx=-yo-zo
    dy=xo+a*yo
    dz=b+zo*(xo-c)

    #Output
    return(c(dx,dy,dz))
  }
  ros.sim<-rk(fun=ros.equ,parameter=parameter,setting=setting,time=time)

  # To avoid an unstable system
  if (nrow(ros.sim) != nrow(ros.sim[complete.cases(ros.sim),])){stop("should simulate it again due to some values of the parameters or initial conditions have leaded to an unstable system that tend to infinity.")}

  # Output
  return(ros.sim[-c(1:(.1*n)),])
}

################################################################################
# We have implemented the classical Runge–Kutta method (RK4) in order to generate
# time-series data from continuous-time dynamical system as the Rössler system.
rk<-function(fun,parameter,setting,time){

  # Checks
  if (!is.function(fun)){stop("'fun' should be a function containing a continuous-time dynamical system")}
  if (!is.vector(parameter)){stop("'parameter' should be a numeric vector denoting the parameter set")}
  if (!is.vector(setting)){stop("'setting' should be a numeric vector denoting the setting set")}
  if (!is.vector(time)){stop("'time' should be a numeric vector denoting the time-lapse and the time-step")}

  # Settings
  n                 <-  length(time)
  out               <-  matrix(ncol=length(setting), nrow=n)
  out[1,]           <-  setting
  h                 <-  time[2]-time[1]

  # Runge-Kutta method (RK4)
  for (i in 2:n){
    k1            <-  h*fun(parameter,out[i-1,],time[i-1])
    k2            <-  h*fun(parameter,out[i-1,]+k1/2,time[i-1]+h/2)
    k3            <-  h*fun(parameter,out[i-1,]+k2/2,time[[i-1]]+h/2)
    k4            <-  h*fun(parameter,out[i-1,]+k3,time[[i-1]]+h)
    out[i,]       <-  out[i-1,]+(1/6)*(k1+2*k2+2*k3+k4)
  }

  # Output
  return(cbind("x"=stats::as.ts(out[,1]),"y"=stats::as.ts(out[,2]),"z"=stats::as.ts(out[,3])))
}






