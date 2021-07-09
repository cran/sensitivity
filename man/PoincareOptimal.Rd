\name{PoincareOptimal}
\alias{PoincareOptimal}

\title{Optimal Poincare constants for Derivative-based Global Sensitivity Measures (DGSM)}

\description{
  A DGSM is a sensitivity index relying on the integral (over the space domain of the input variables) 
  of the squared derivatives of a model output with respect to one model input variable. 
  The product between a DGSM and a Poincare Constant (Roustant et al., 2014: Roustant et al., 2017),
  on the type of probability distribution of the input variable, gives an   upper bound of the total 
  Sobol' index corresponding to the same input (Lamboni et al., 2013; Kucherenko and Iooss, 2016).
  
  This function provides the optimal Poincare constant as explained in Roustant et al. (2017).
  It solves numerically the spectral problem corresponding to the Poincare inequality, with 
  Neumann conditions. The differential equation is  f'' - V'f'= - lambda f with f'(a) = f'(b) = 0.
  In addition, all the spectral decomposition can be returned by the function. 
  The eigenvalues are sorted in ascending order, starting from zero. 
  The information corresponding to the optimal constant is thus given in the second column.
  
  IMPORTANT: This program is useless for the two following input variable distributions:
  \itemize{
    \item uniform on \eqn{[min,max]}{[min,max]} interval: The optimal Poincare constant is \eqn{\frac{(max-min)^2}{pi^2}}{(max-min)^2/(pi^2)}.
    \item normal with a standard deviation \eqn{sd}{sd}: The optimal Poincare constant is \eqn{sd^2}{sd^2}.
    }
}

\usage{
PoincareOptimal(distr=list("unif",c(0,1)), min=NULL, max=NULL, 
                n = 500, method = c("quadrature", "integral"), only.values = TRUE, 
                der = FALSE, plot = FALSE, ...)
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
  \item{der}{if TRUE, compute the eigenfunction derivatives (default value is FALSE)}
  \item{plot}{logical:if TRUE and only.values=FALSE, plots a minimizer of the Rayleigh ratio (default value is FALSE)}
  \item{...}{additional arguments}
}

\value{
  \code{PoincareOptimal} returns a list containing:
  \item{opt}{the optimal Poincare constant}
  \item{values}{the eigenvalues in increasing order, starting from 0. Thus, the second one is the spectral gap, equal to the inverse of the Poincare constant}
  \item{vectors}{the values of eigenfunctions at \code{knots}}
  \item{der}{the values of eigenfunction derivatives at \code{knots}}
  \item{knots}{a sequence of length \code{n} formed by equally spaced real numbers in the support of the probability distribution, used for discretization}
}

\details{
For the uniform, normal, triangular and Gumbel distributions, the optimal constants are computed on the standardized correponding distributions (for a better numerical efficiency). In these cases, the return optimal constant and eigenvalues correspond to original distributions.
}

\references{
O. Roustant, F. Barthe and B. Iooss, 
Poincare inequalities on intervals - application to sensitivity analysis,
Electronic Journal of Statistics, Vol. 11, No. 2, 3081-3119, 2017.

O Roustant, F. Gamboa, B Iooss. Parseval inequalities and lower bounds 
# for variance-based sensitivity indices. 2019. hal-02140127
}

\author{
Olivier Roustant and Bertrand Iooss
}

\seealso{
  \code{\link{PoincareConstant}, \link{PoincareChaosSqCoef}}
}

\examples{


# uniform on [a, b]
a <- -1 ; b <- 1
out <- PoincareOptimal(distr = list("unif", a, b))
cat("Poincare constant (theory -- estimated):", (b-a)^2/pi^2, "--", out$opt, "\n")

# truncated standard normal on [-1, 1]
# the optimal Poincare constant is then equal to 1/3,
# as -1 and 1 are consecutive roots of the 2nd Hermite polynomial X*X - 1.
out <- PoincareOptimal(distr = dnorm, min = -1, max = 1, 
                       plot = TRUE, only.values = FALSE)
cat("Poincare constant (theory -- estimated):", 1/3, "--", out$opt, "\n")

\donttest{
# truncated standard normal on [-1.87, +infty]
out <- PoincareOptimal(distr = list("norm", 0, 1), min = -1.87, max = 5, 
                       method = "integral", n = 500)
print(out$opt)

# truncated Gumbel(0,1) on [-0.92, 3.56]
library(evd)
out <- PoincareOptimal(distr = list("gumbel", 0, 1), min = -0.92, max = 3.56, 
                       method = "integral", n = 500)
print(out$opt)

# symetric triangular [-1,1]
library(triangle)
out <- PoincareOptimal(distr = list("triangle", -1, 1, 0), min = NULL, max = NULL)
cat("Poincare constant (theory -- estimated):", 0.1729, "--", out$opt, "\n")


# Lognormal distribution
out <- PoincareOptimal(distr = list("lognorm", 1, 2), min = 3, max = 10, 
                       only.values = FALSE, plot = TRUE, method = "integral")
print(out$opt)


## -------------------------------

## Illustration for eigenfunctions on the uniform distribution
## (corresponds to Fourier series)
b <- 1
a <- -b
out <- PoincareOptimal(distr = list("unif", a, b), 
                       only.values = FALSE, der = TRUE, method = "quad")

# Illustration for 3 eigenvalues

par(mfrow = c(3,2))
eigenNumber <- 1:3 # eigenvalue number
for (k in eigenNumber[1:3]){ # keep the 3 first ones (for graphics)
  plot(out$knots, out$vectors[, k + 1], type = "l", 
       ylab = "", main = paste("Eigenfunction", k), 
       xlab = paste("Eigenvalue:", round(out$values[k+1], digits = 3)))
  sgn <- sign(out$vectors[1, k + 1])
  lines(out$knots, sgn * sqrt(2) * cos(pi * k * (out$knots/(b-a) + 0.5)), 
        col = "red", lty = "dotted")
  
  plot(out$knots, out$der[, k + 1], type = "l", 
       ylab = "", main = paste("Eigenfunction derivative", k), 
       xlab = "")
  sgn <- sign(out$vectors[1, k + 1])
  lines(out$knots, - sgn * sqrt(2) / (b-a) * pi * k * sin(pi * k * (out$knots/(b-a) + 0.5)), 
        col = "red", lty = "dotted")
}


# how to create a function for one eigenfunction and eigenvalue,
# given N values 
eigenFun <- approxfun(x = out$knots, y = out$vectors[, 2])
eigenDerFun <- approxfun(x = out$knots, y = out$der[, 2])
x <- runif(n = 3, min = -1/2, max = 1/2)
eigenFun(x)
eigenDerFun(x)

}
}
