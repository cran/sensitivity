\name{sobolmara}
\alias{sobolmara}
\alias{tell.sobolmara}
\alias{print.sobolmara}
\alias{plot.sobolmara}

\title{Monte Carlo Estimation of Sobol' Indices via matrix permutations}

\description{ \code{sobolmara} implements the Monte Carlo estimation of
   the first-order Sobol' sensitivity indices using the formula of Mara and Joseph (2008), called the
   Mara estimator.
   This method allows the estimation of all first-order p indices at a cost of 
   2N model calls (the random sample size), then independently of p (the number of inputs). }

\usage{
sobolmara(model = NULL, X1, \dots)
\method{tell}{sobolmara}(x, y = NULL, return.var = NULL, \dots)
\method{print}{sobolmara}(x, \dots)
\method{plot}{sobolmara}(x, ylim = c(0, 1), \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{X1}{the random sample.}
  \item{x}{a list of class \code{"sobolEff"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{return.var}{a vector of character strings giving further
    internal variables names to store in  the output object \code{x}.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
}

\value{
  \code{sobolmara} returns a list of class \code{"sobolmara"}, containing all
  the input arguments detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments.}
  \item{y}{a vector of model responses.}
  \item{S}{the estimations of the Sobol' sensitivity indices.}

}

\details{
	The estimator used by sobolmara is based on rearragement of a unique matrix via random permutations (see Mara and Joseph, 2008). Bootstrap confidence intervals are not available.
}

\references{
	Mara, T. and Joseph, O.R. (2008), \emph{Comparison of some efficient methods to evaluate the main effect of computer model factors}, Journal of Statistical Computation and Simulation, 78:167--178
}

\author{
Bertrand Iooss
}

\seealso{
  \code{\link{sobolroalhs}, \link{sobol}, \link{sobol2002}, \link{sobol2007}, \link{soboljansen}, \link{sobolmartinez}, \link{sobolEff},\link{sobolMultOut}}
}

\examples{
# Test case : the non-monotonic Sobol g-function

# The method of sobolmara requires 1 sample
# (there are 8 factors, all following the uniform distribution on [0,1])
n <- 1000
X1 <- data.frame(matrix(runif(8 * n), nrow = n))

# sensitivity analysis
x <- sobolmara(model = sobol.fun, X1 = X1)
print(x)
plot(x)
}

\keyword{design}
