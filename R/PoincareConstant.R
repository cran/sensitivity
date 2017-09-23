# Poincare Constant computation for Derivative-based Global Sensitivity Measures (DGSM)
# using log-concave case formula, double exponential transport or logistic transport
#
# Authors: Jana Fruth (2014), Bertrand Iooss and Olivier Roustant (2016)
#
# References: 
#
# O. Roustant, J. Fruth, B. Iooss and S. Kuhnt,
# Crossed-derivative-based sensitivity measures for interaction screening, 
# Mathematics and Computers in Simulation, 105:105-118, 2014
#
#O. Roustant, F. Barthe and B. Iooss,
#Poincare inequalities on intervals - application to sensitivity analysis, 
#Electronic Journal of Statistics, Vol. 11, No. 2, 3081-3119, 2017

PoincareConstant <- function(dfct=dnorm, qfct=qnorm, pfct=pnorm, 
                             logconcave=FALSE, transport="logistic", optimize.interval=c(-100, 100), 
                             truncated=FALSE, min=0, max=1, ...){
  
  if (logconcave == TRUE){
    if (transport == "logistic") warning("Double-exponential transport has been used here instead of the logistic one, since the analytical formula based on the log-concave law cases is a subcase of the double-exponential transport.")

    if (truncated == FALSE) res <- 1/dfct(qfct(0.5,...),...)^2
    if (truncated == TRUE){
      res <- (pfct(max,...)- pfct(min,...))^2 /
        (dfct(qfct((pfct(min,...)+pfct(max,...))/2,...),...))^2
    }
  }
  if (logconcave == FALSE){
    if (transport == "double_exp"){
      fct <- function(x){
        if (truncated == FALSE){
          cdf.at.x <- pfct(x, ...)
          density.at.x <- dfct(x, ...)
        }
        if (truncated == TRUE){
          cdf.at.x <- pfct(x, min=min, max=max, ...)
          density.at.x <- dfct(x, min=min, max=max, ...)
        }
        apply(cbind(cdf.at.x, 1-cdf.at.x),1,min)/(density.at.x)
      }
    }
    if (transport == "logistic"){
      fct <- function(x){
        if (truncated == FALSE){
          cdf.at.x <- pfct(x, ...)
          density.at.x <- dfct(x, ...)
        }
        if (truncated == TRUE){
          cdf.at.x <- pfct(x, min=min, max=max, ...)
          density.at.x <- dfct(x, min=min, max=max, ...)
        }
        (cdf.at.x) * (1-cdf.at.x)/(density.at.x)
      }
    }
    
    c1 <- optimize(f=fct, interval=optimize.interval, maximum=TRUE)$objective
    res <- 4*c1^2
  }
  print(res)
}
