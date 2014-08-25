\name{PoincareConstant}
\alias{PoincareConstant}

\title{Poincare constants for Derivative-based Global Sensitivity Measures (DGSM)}

\description{
  A DGSM is the product between a Poincare Constant (Roustant et al., 2014) and the integral 
  (over the space domain of the input variables) of the squared derivatives of the model output 
  with respect to one input variable. The DGSM is a maximal bound of the total Sobol' index 
  corresponding to the same input (Lamboni et al., 2013).
  
  This DGSM depends on the type of probability distribution of the input variable.
  In the particular case of log-concave distribution, analytical formulas are available by
  the way of the median value (Lamboni et al., 2013).
  For truncated log-concave distributions, different formulas are available (Roustant et al., 2014).
  For general non-truncated distributions (including the non log-concave case), the Poincare constant
  is computed via a relatively simple optimization process (Lamboni et al., 2013).
  
  IMPORTANT: This program is useless for the two following input variable distributions:
  \itemize{
    \item uniform on \eqn{[min,max]}{[min,max]} interval: The optimal Poincare constant is \eqn{\frac{(max-min)^2}{pi^2}}{(max-min)^2/(pi^2)}.
    \item normal with a standard deviation \eqn{sd}{sd}: The optimal Poincare constant is \eqn{sd^2}{sd^2}.
  }
}

\usage{
PoincareConstant(densityfct=dnorm, qfct=qnorm, cdfct,
                 truncated=FALSE, min=0, max=1, 
                 logconcave=TRUE, optimize.interval=c(-100, 100), ...)
}

\arguments{
  \item{densityfct}{the probability density function of the input variable}
  \item{qfct}{the quantile function of the input variable}
  \item{cdfct}{the distribution function of the input variable}
  \item{truncated}{logical value: TRUE for a truncated distribution. Default value is FALSE}
  \item{min}{the minimal bound in the case of a truncated distribution}
  \item{max}{the maximal bound in the case of a truncated distribution}
  \item{logconcave}{logical value: TRUE (default value) for a log-concave distribution}
  \item{optimize.interval}{In the non-log concave case, a vector containing the end-points of the interval to be searched for the maximum of the function to be optimized}
  \item{...}{additional arguments}
}

\value{
  \code{PoincareConstant} returns the value of the Poincare constant.
}

\references{
O. Roustant, J. Fruth, B. Iooss and S. Kuhnt,
Crossed-derivative-based sensitivity measures for interaction screening, 
Mathematics and Computers in Simulation, 105:105-118, 2014.

M. Lamboni, B. Iooss, A-L. Popelin and F. Gamboa,
Derivative-based global sensitivity measures: General links with Sobol' indices and 
numerical tests, Mathematics and Computers in Simulation, 87:45-54, 2013.
}

\author{
Jana Fruth and Bertrand Iooss
}

\examples{

# Exponential law (log-concave)
PoincareConstant(dexp,qexp,rate=1)

# Weibull law (non log-concave)
PoincareConstant(dweibull,cdfct=pweibull, logconcave=FALSE, 
optimize.interval=c(0, 15), shape=1, scale=1)

\dontrun{
# Triangular law (log-concave)
library(triangle)
PoincareConstant(dtriangle, qtriangle, a=49, b=51, c=50)

# Truncated Gumbel law (log-concave)
library(evd)
PoincareConstant(dgumbel, qgumbel, pgumbel, truncated=TRUE, 
min=500, max=3000, loc=1013.0, scale=558.0) 
}

}
