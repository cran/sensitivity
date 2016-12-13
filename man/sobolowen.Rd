\name{sobolowen}
\alias{sobolowen}
\alias{tell.sobolowen}
\alias{print.sobolowen}
\alias{plot.sobolowen}

\title{Monte Carlo Estimation of Sobol' Indices (improved formulas of Owen (2013)}

\description{
  \code{sobolowen} implements the Monte Carlo estimation of
  the Sobol' indices for both first-order and total indices at the same
  time (alltogether \eqn{2p}{2p} indices). Take as input 3 independent matrices.
  These are called the Owen estimators.
}

\usage{
sobolowen(model = NULL, X1, X2, X3, nboot = 0, conf = 0.95, varest = 2, \dots)
\method{tell}{sobolowen}(x, y = NULL, return.var = NULL, varest = 2, \dots)
\method{print}{sobolowen}(x, \dots)
\method{plot}{sobolowen}(x, ylim = c(0, 1), \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{X1}{the first random sample.}
  \item{X2}{the second random sample.}
  \item{X3}{the third random sample.}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level for bootstrap confidence intervals.}
  \item{varest}{choice for the variance estimator for the denominator of
  the Sobol' indices. varest=1 is for a classical estimator.
  varest=2 (default) is for the estimator proposed in Janon et al. (2012).}
  \item{x}{a list of class \code{"sobol"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{return.var}{a vector of character strings giving further
    internal variables names to store in  the output object \code{x}.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called}
}

\value{
  \code{sobolowen} returns a list of class \code{"sobolowen"}, containing all
  the input arguments detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments.}
  \item{y}{the response used}
  \item{V}{the estimations of Variances of the Conditional Expectations
    (VCE) with respect to each factor and also with respect to the
  complementary set of each factor ("all but \eqn{X_i}{Xi}").}
  \item{S}{the estimations of the Sobol' first-order indices.}
  \item{T}{the estimations of the Sobol' total sensitivity indices.}
  
  Users can ask more ouput variables with the argument
  \code{return.var} (for example, bootstrap outputs \code{V.boot},
  \code{S.boot} and \code{T.boot}).
}

\references{
A. Owen, 2013, \emph{Better estimations of small Sobol' sensitivity indices}, 
ACM Transactions on Modeling and Computer Simulations (TOMACS), 23(2), 11.

Janon, A., Klein T., Lagnoux A., Nodet M., Prieur C. (2012), Asymptotic 
normality and efficiency of two Sobol index estimators. Accepted in 
ESAIM: Probability and Statistics.
}

\author{
Taieb Touati and Bernardo Ramos
}

\seealso{
  \code{\link{sobol}, \link{sobol2002}, \link{sobolSalt}, \link{sobol2007}, \link{soboljansen}, \link{sobolmartinez}, \link{sobolEff}, \link{sobolmara}, \link{sobolGP}}
}

\examples{
# Test case : the non-monotonic Sobol g-function

# The method of sobolowen requires 3 samples
# There are 8 factors, all following the uniform distribution
# on [0,1]

library(boot)
n <- 1000
X1 <- data.frame(matrix(runif(8 * n), nrow = n))
X2 <- data.frame(matrix(runif(8 * n), nrow = n))
X3 <- data.frame(matrix(runif(8 * n), nrow = n))

# sensitivity analysis

\dontrun{
x <- sobolowen(model = sobol.fun, X1, X2, X3, nboot = 100)
print(x)
plot(x)
}
}

\keyword{design}
