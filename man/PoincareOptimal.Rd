\name{PoincareOptimal}
\alias{PoincareOptimal}

\title{Optimal Poincare constants for Derivative-based Global Sensitivity Measures (DGSM)}

\description{
  A DGSM is a sensitivity index relying on the integral (over the space domain of the input variables) 
  of the squared derivatives of a model output with respect to one model input variable. 
  The product between a DGSM and a Poincare Constant (Roustant et al., 2014: Roustant et al., 2016),
  on the type of probability distribution of the input variable, gives an   upper bound of the total 
  Sobol' index corresponding to the same input (Lamboni et al., 2013; Kucherenko and Iooss, 2016).
  
  This function provides the optimal Poincare constant as explained in Roustant et al. (2016).
  It solves numerically the spectral problem corresponding to the Poincare inequality, with 
  Neumann conditions. The differential equation is  f'' - V'f'= - lambda f with f'(a) = f'(b) = 0.
  In addition, all the spectral decomposition can be returned by the function. The information
  corresponding to the optimal constant is given in the second to last column.
  
  IMPORTANT: This program is useless for the two following input variable distributions:
  \itemize{
    \item uniform on \eqn{[min,max]}{[min,max]} interval: The optimal Poincare constant is \eqn{\frac{(max-min)^2}{pi^2}}{(max-min)^2/(pi^2)}.
    \item normal with a standard deviation \eqn{sd}{sd}: The optimal Poincare constant is \eqn{sd^2}{sd^2}.
    }
}

\usage{
PoincareOptimal(distr=list("unif",c(0,1)), min=NULL, max=NULL, n = 500, 
    method = c("quadrature", "integral"), only.values = TRUE, plot = FALSE, ...)
}

\arguments{
  \item{distr}{a list or a function corresponding to the probability distribution. 
  \itemize{
  \item If it is a list, it contains the name of the R distribution of the variable and its parameters. Possible choices are: "unif" (uniform), "norm" (normal), "exp" (exponential), "triangle" (triangular from package triangle), "gumbel" (from package evd), "beta", "gamma", "weibull" and "lognorm" (lognormal). The values of the distribution parameters have to be passed in arguments in the same order than the corresponding R function.
  \item If it is a function, it corresponds to the pdf. Notice that the normalizing constant has no impact on the computation of the optimal Poincare constant and can be ommitted.
  }}
  \item{min}{see below}
  \item{max}{[min,max]: interval on which the distribution is truncated. Choose low and high quantiles in case of unbounded distribution. Choose NULL for uniform and triangular distributions}
  \item{n}{number of discretization steps}
  \item{method}{method of integration: "quadrature" (default value) uses the trapez quadrature (close and quicker), "integral" is longer but does not make any approximation}
  \item{only.values}{if TRUE, only eigen values are computed and returned, otherwise both eigenvalues and   eigenvectors are returned (default value is TRUE)}
  \item{plot}{logical:if TRUE and only.values=FALSE, plots a minimizer of the Rayleigh ratio (default value is FALSE)}
  \item{...}{additional arguments}
}

\value{
  \code{PoincareOptimal} returns a list containing:
  \item{opt}{the optimal Poincare constant}
  \item{values}{the eigen values}
  \item{vectors}{the eigen vectors}
}

\details{
For the uniform, normal, triangular and Gumbel distributions, the optimal constants are computed on the standardized correponding distributions (for a better numerical efficiency). In these cases, the return optimal constant and eigen values correspond to original distributions, while the eigen vectors are not rescaled.
}

\references{
O. Roustant, F. Barthe and B. Iooss, 
Poincare inequalities on intervals - application to sensitivity analysis,
Submitted, 2016. 
}

\author{
Olivier Roustant and Bertrand Iooss
}

\seealso{
  \code{\link{PoincareConstant}}
}

\examples{


# uniform [0,1]
out <- PoincareOptimal(distr=list("unif",0,1))
print(out$opt)

# truncated standard normal on [-1, 1]
out <- PoincareOptimal(distr=dnorm, min=-1, max=1, plot=TRUE, only.values=FALSE)
print(out$opt)

\dontrun{
# truncated standard normal on [-1.87, +infty]
out <- PoincareOptimal(distr=list("norm",0,1), min=-1.87, max=5, method="integral", n=500)
print(out$opt)

# truncated Gumbel(0,1) on [-0.92, 3.56]
out <- PoincareOptimal(distr=list("gumbel",0,1), min=-0.92, max=3.56, method="integral", n=500)
print(out$opt)

# symetric triangular [-1,1]
out <- PoincareOptimal(distr=list("triangle",-1,1,0), min=NULL, max=NULL)
print(out$opt)

# Lognormal distribution
out <- PoincareOptimal(distr=list("lognorm",1,2), min=3, max=10, only.values=FALSE,plot=TRUE, 
  method="integral")
print(out$opt)

}

}
