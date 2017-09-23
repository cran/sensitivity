\name{sobolMultOut}
\alias{sobolMultOut}
\alias{print.sobolMultOut}
\alias{plot.sobolMultOut}

\title{Monte Carlo Estimation of Aggregated Sobol' Indices for multiple and functional outputs}

\description{
  \code{sobolMultOut} implements the aggregated Sobol' indices for 
  multiple outputs. It consists in averaging all the Sobol indices weighted
  by the variance of their corresponding output. Moreover, this function computes and plots
  the functional (unidimensional) Sobol' indices for functional (unidimensional)
  model output. Sobol' indices for both first-order and total indices are estimated 
  by Monte Carlo formulas.
}

\usage{
sobolMultOut(model = NULL, q = 1, X1, X2, MCmethod = "sobol", 
                         plotFct=FALSE, \dots)
\method{print}{sobolMultOut}(x, \dots)
\method{plot}{sobolMultOut}(x, ylim = c(0, 1), \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{q}{dimension of the model output vector.}
  \item{X1}{the first random sample.}
  \item{X2}{the second random sample.}
  \item{MCmethod}{a character string specifying the Monte-Carlo procedure used 
  to estimate the Sobol indices. The avaible methods are : \code{"sobol"}, 
  \code{"sobol2002"}, \code{"sobol2007"}, \code{"soboljansen"}, \code{sobolmara}
  and \code{sobolGP}.}
  \item{plotFct}{if TRUE, 1D functional Sobol indices are computed and plotted 
  in an external window (default=FALSE).}
  \item{x}{a list of class \code{MCmethod} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{ylim}{y-coordinate plotting limits.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called}
}

\value{
  \code{sobolMultOut} returns a list of class \code{MCmethod}, containing all
  its input arguments, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments.}
  \item{y}{the response used}
  \item{V}{the estimations of the aggregated Variances of the Conditional Expectations
    (VCE) with respect to each factor and also with respect to the
  complementary set of each factor ("all but \eqn{X_i}{Xi}").}
  \item{S}{the estimations of the aggregated Sobol' first-order indices.}
  \item{T}{the estimations of the aggregated Sobol' total sensitivity indices.}
  \item{Sfct}{the estimations of the functional Sobol' first-order indices 
  (if PlotFct=TRUE).}
  \item{Tfct}{the estimations of the functional Sobol' total sensitivity indices
  (if PlotFct=TRUE).}
}

\details{
For this function, there are several gaps: the bootstrap estimation of confidence 
intervals is not avalaible and the tell function does not work.
}

\references{
M. Lamboni, H. Monod and D. Makowski, 2011, \emph{Multivariate sensitivity analysis to measure global contribution of input factors in dynamic models}, Reliability Engineering and System Safety, 96:450-459.

F. Gamboa, A. Janon, T. Klein and A. Lagnoux, 2014, \emph{Sensitivity indices for multivariate outputs}, Electronic Journal of Statistics, 8:575-603.
}

\author{
Bertrand Iooss
}

\seealso{
  \code{\link{sobol}, \link{sobol2002}, \link{sobol2007}, \link{soboljansen},
  \link{sobolmara}, \link{sobolGP}}
}

\examples{
\dontrun{
# Functional toy function: Arctangent temporal function (Auder, 2011)
# X: input matrix (in [-7,7]^2)
# q: number of discretization steps of [0,2pi] interval
# output: vector of q values

atantemp <- function(X, q = 100){
  
  n <- dim(X)[[1]]
  t <- (0:(q-1)) * (2*pi) / (q-1)
  
  res <- matrix(0,ncol=q,nrow=n)
  for (i in 1:n) res[i,] <- atan(X[i,1]) * cos(t) + atan(X[i,2]) * sin(t)

  return(res)  
}

# Tests functional toy fct 

y0 <- atantemp(matrix(c(-7,0,7,-7,0,7),ncol=2))
#plot(y0[1,],type="l")
#apply(y0,1,lines)

n <- 100
X <- matrix(c(runif(2*n,-7,7)),ncol=2)
y <- atantemp(X)
x11()
plot(y0[2,],ylim=c(-2,2),type="l")
apply(y,1,lines)

# Sobol indices computations

n <- 1000
X1 <- data.frame(matrix(runif(2*n,-7,7), nrow = n))
X2 <- data.frame(matrix(runif(2*n,-7,7), nrow = n))

x11()
sa <- sobolMultOut(model=atantemp, q=100, X1, X2, 
                   MCmethod="soboljansen", plotFct=T)
print(sa)
x11()
plot(sa)

}
}

\keyword{design}
