\name{sobolEff}
\alias{sobolEff}
\alias{tell.sobolEff}
\alias{print.sobolEff}
\alias{plot.sobolEff}

\title{Monte Carlo Estimation of Sobol' Indices (formulas of Janon-Monod)}

\description{ \code{sobolEff} implements the Monte Carlo estimation of the Sobol' sensitivity indices using the asymptotically efficient formulas in section 4.2.4.2 of Monod et al. (2006). Either all first-order indices or all total-effect indices are estimated at a cost of \eqn{N \times (p+1)}{N*(p+1)} model calls or all closed second-order indices are estimated at a cost of \eqn{N \times p \choose 2)}{N*(p*(p-1)/2+1)} model calls.}

\usage{
sobolEff(model = NULL, X1, X2, order=1, nboot = 0, conf = 0.95, \dots)
\method{tell}{sobolEff}(x, y = NULL, \dots)
\method{print}{sobolEff}(x, \dots)
\method{plot}{sobolEff}(x, ylim = c(0, 1), \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{X1}{the first random sample.}
  \item{X2}{the second random sample.}
  \item{order}{an integer specifying the indices to estimate: 0 for total effect indices,1 for first-order indices and 2 for closed second-order indices.}
  \item{nboot}{the number of bootstrap replicates, or zero to use asymptotic 
  standard deviation estimates given in Janon et al. (2012).}
  \item{conf}{the confidence level for confidence intervals.}
  \item{x}{a list of class \code{"sobolEff"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
}

\value{
  \code{sobolEff} returns a list of class \code{"sobolEff"}, containing all
  the input arguments detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments.}
  \item{y}{a vector of model responses.}
  \item{S}{the estimations of the Sobol' sensitivity indices.}

}

\details{
		 The estimator used by sobolEff is defined in Monod et al. (2006), 
     Section 4.2.4.2 and studied under the name T_N in Janon et al. (2012). 
     This estimator is good for large first-order indices.
}

\references{
		 Monod, H., Naud, C., Makowski, D. (2006), Uncertainty and sensitivity 
     analysis for crop models in Working with Dynamic Crop Models: Evaluation, 
     Analysis, Parameterization, and Applications, Elsevier.

A. Janon, T. Klein, A. Lagnoux, M. Nodet, C. Prieur (2014), \emph{Asymptotic normality and efficiency of two Sobol index estimators}, ESAIM: Probability and Statistics, 18:342-364.

}

\author{
Alexandre Janon, Laurent Gilquin
}

\seealso{
  \code{\link{sobol}, \link{sobol2002}, \link{sobolSalt}, \link{sobol2007}, \link{soboljansen}, \link{sobolmartinez},
  \link{sobolSmthSpl}}
}

\examples{
# Test case : the non-monotonic Sobol g-function

# The method of sobol requires 2 samples
# (there are 8 factors, all following the uniform distribution on [0,1])
n <- 1000
X1 <- data.frame(matrix(runif(8 * n), nrow = n))
X2 <- data.frame(matrix(runif(8 * n), nrow = n))

# sensitivity analysis
x <- sobolEff(model = sobol.fun, X1 = X1, X2 = X2, nboot = 0)
print(x)
}

\keyword{design}