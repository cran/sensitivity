\name{PoincareConstant}
\alias{PoincareConstant}

\title{Poincare constants for Derivative-based Global Sensitivity Measures (DGSM)}

\description{
  A DGSM is a sensitivity index relying on the integral (over the space domain of the input variables) 
  of the squared derivatives of a model output with respect to one model input variable. 
  The product between a DGSM and a Poincare Constant (Roustant et al., 2014: Roustant et al., 2016) 
  gives an upper bound of the total Sobol' index corresponding to the same input 
  (Lamboni et al., 2013; Kucherenko and Iooss, 2016).
  
  This Poincare constant depends on the type of probability distribution of the input variable.
  In the particular case of log-concave distribution, analytical formulas are available for
  double-exponential transport by the way of the median value (Lamboni et al., 2013). For truncated 
  log-concave distributions, different formulas are available (Roustant et al., 2014). For general 
  distributions (truncated or not), some Poincare constants can be computed via a relatively simple
  optimization process using different formula coming from transport inequalities (Roustant et al., 2016).
  
  Notice that the analytical formula based on the log-concave law cases is a subcase of the 
  double-exponential transport. In all cases, with this function, the smallest constant is obtained using
  the logistic transport formula. \code{\link{PoincareOptimal}} allows to obtained the best (optimal)
  constant using another (spectral) method.
  
  IMPORTANT: This program is useless for the two following input variable distributions:
  \itemize{
    \item uniform on \eqn{[min,max]}{[min,max]} interval: The optimal Poincare constant is \eqn{\frac{(max-min)^2}{pi^2}}{(max-min)^2/(pi^2)}.
    \item normal with a standard deviation \eqn{sd}{sd}: The optimal Poincare constant is \eqn{sd^2}{sd^2}.
  }
}

\usage{
PoincareConstant(dfct=dnorm, qfct=qnorm, pfct=pnorm,
                 logconcave=FALSE, transport="logistic", optimize.interval=c(-100, 100),
                 truncated=FALSE, min=0, max=1,  ...)
}

\arguments{
  \item{dfct}{the probability density function of the input variable}
  \item{qfct}{the quantile function of the input variable}
  \item{pfct}{the distribution function of the input variable}
  \item{logconcave}{logical value: TRUE for a log-concave distribution (analyical formula will be used). Requires argument 'dfct' and 'qfct'. FALSE (default value) means that the calculations will be performed using transport-based formulas (applicable for log-concave and non-log concave cases)}
  \item{transport}{If logconcave=FALSE, choice of the transport inequalities to be used: "double_exp" (default value) for double exponential transport and "logistic" for logistic transport". Requires argument 'dfct' and 'pfct'}
  \item{optimize.interval}{In the transport-based case (logconcave=FALSE), a vector containing the end-points of the interval to be searched for the maximum of the function to be optimized}
  \item{truncated}{logical value: TRUE for a truncated distribution. Default value is FALSE}
  \item{min}{the minimal bound in the case of a truncated distribution}
  \item{max}{the maximal bound in the case of a truncated distribution}
  \item{...}{additional arguments}
}

\value{
  \code{PoincareConstant} returns the value of the Poincare constant.
}

\details{
	In the case of truncated distributions (truncated=TRUE), in addition to the min and max arguments: - the truncated distribution name has to be passed in the 'dfct' and 'pfct' arguments if logconcave=FALSE, - the non-truncated distribution name has to be passed in the 'dfct' and 'qfct' arguments if logconcave=TRUE. Moreover, if min and max are finite, optimize.interval is required to be defined as c(min,max).
}

\references{
S. Kucherenko and B. Iooss,
Derivative-based global sensitivity measures,
In: R. Ghanem, D. Higdon and H. Owhadi (eds.), Handbook of Uncertainty Quantification, 
2016.

M. Lamboni, B. Iooss, A-L. Popelin and F. Gamboa,
Derivative-based global sensitivity measures: General links with Sobol' indices and 
numerical tests, Mathematics and Computers in Simulation, 87:45-54, 2013.

O. Roustant, F. Barthe and B. Iooss,
Poincare inequalities on intervals - application to sensitivity analysis, 
Submitted, 2016. 

O. Roustant, J. Fruth, B. Iooss and S. Kuhnt,
Crossed-derivative-based sensitivity measures for interaction screening, 
Mathematics and Computers in Simulation, 105:105-118, 2014.
}

\author{
Jana Fruth, Bertrand Iooss and Olivier Roustant
}

\seealso{
  \code{\link{PoincareOptimal}}
}

\examples{

# Exponential law (log-concave)
PoincareConstant(dfct=dexp,qfct=qexp,pfct=NULL,rate=1,logconcave=TRUE) # log-concave assumption
PoincareConstant(dfct=dexp,qfct=NULL,pfct=pexp,rate=1,optimize.interval=c(0, 15)) 
          # logistic transport approach

# Weibull law (log-concave)
PoincareConstant(dfct=dweibull,qfct=NULL,pfct=pweibull,optimize.interval=c(0, 15),shape=1,scale=1) 
          # logistic transport approach

\dontrun{
# Triangular law (log-concave)
library(triangle)
PoincareConstant(dfct=dtriangle, qfct=qtriangle, pfct=NULL, a=-1, b=1, c=0, logconcave=TRUE) 
          # log-concave assumption
PoincareConstant(dfct=dtriangle, qfct=NULL, pfct=ptriangle, a=-1, b=1, c=0, 
  transport="double_exp", optimize.interval=c(-1,1)) # Double-exponential transport approach
PoincareConstant(dfct=dtriangle, qfct=NULL, pfct=ptriangle, a=-1, b=1, c=0, 
  optimize.interval=c(-1,1)) # Logistic transport calculation

# Normal N(0,1) law truncated on [-1.87,+infty]
PoincareConstant(dfct=dnorm, qfct=qnorm, pfct=pnorm, mean=0, sd=1, logconcave=TRUE, 
  transport="double_exp", truncated=TRUE, min=-1.87, max=999) # log-concave assumption
PoincareConstant(dfct=dtnorm, qfct=qtnorm, pfct=ptnorm, mean=0, sd=1, truncated=TRUE, 
  min=-1.87, max=999,   transport="double_exp", optimize.interval=c(-1.87,20)) 
            # Double-exponential transport approach
PoincareConstant(dfct=dtnorm, qfct=qtnorm, pfct=ptnorm, mean=0, sd=1, truncated=TRUE, 
  min=-1.87, max=999, optimize.interval=c(-1.87,20)) # Logistic transport approach


# Gumbel law (log-concave)
library(evd)
PoincareConstant(dfct=dgumbel, qfct=qgumbel, pfct=NULL, loc=0, scale=1, logconcave=TRUE,
  transport="double_exp") # log-concave assumption
PoincareConstant(dfct=dgumbel, qfct=NULL, pfct=pgumbel, loc=0, scale=1, 
  transport="double_exp", optimize.interval=c(-3,20)) # Double-exponential transport approach
PoincareConstant(dfct=dgumbel, qfct=qgumbel, pfct=pgumbel, loc=0, scale=1, 
  optimize.interval=c(-3,20)) # Logistic transport approach

# Truncated Gumbel law (log-concave)
PoincareConstant(dfct=dgumbel, qfct=qgumbel, pfct=pgumbel, loc=0, scale=1, logconcave=TRUE, 
  transport="double_exp", truncated=TRUE, min=-0.92, max=3.56) # log-concave assumption
PoincareConstant(dfct=dtgumbel, qfct=NULL, pfct=ptgumbel, loc=0, scale=1, truncated=TRUE, 
  min=-0.92, max=3.56, transport="double_exp", optimize.interval=c(-0.92,3.56)) 
            # Double-exponential transport approach
PoincareConstant(dfct=dtgumbel, qfct=qtgumbel, pfct=ptgumbel, loc=0, scale=1, truncated=TRUE, 
  min=-0.92, max=3.56, optimize.interval=c(-0.92,3.56)) # Logistic transport approach
  
}

}
