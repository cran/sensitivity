\name{sobol2007}
\alias{sobol2007}
\alias{tell.sobol2007}
\alias{print.sobol2007}
\alias{plot.sobol2007}

\title{Monte Carlo Estimation of Sobol' Indices (improved formulas of Sobol et al. (2007) and Saltelli et al. (2010))}

\description{
  \code{sobol2007} implements the Monte Carlo estimation of
  the Sobol' indices for both first-order and total indices at the same
  time (alltogether \eqn{2p}{2p} indices), at a total cost of \eqn{(p+2)
    \times n}{(p + 2) * n} model evaluations. These are called the Mauntz estimators.
}

\usage{
sobol2007(model = NULL, X1, X2, nboot = 0, conf = 0.95, \dots)
\method{tell}{sobol2007}(x, y = NULL, return.var = NULL, \dots)
\method{print}{sobol2007}(x, \dots)
\method{plot}{sobol2007}(x, ylim = c(0, 1), \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{X1}{the first random sample.}
  \item{X2}{the second random sample.}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level for bootstrap confidence intervals.}
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
  \code{sobol2007} returns a list of class \code{"sobol2007"}, containing all
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

\details{
This estimator is good for small first-order and total indices.
}

\references{
I.M. Sobol, S. Tarantola, D. Gatelli, S.S. Kucherenko and W. Mauntz, 2007, \emph{Estimating
the approximation errors when fixing unessential factors in global sensitivity analysis},
Reliability Engineering and System Safety, 92, 957--960.

A. Saltelli, P. Annoni, I. Azzini, F. Campolongo, M. Ratto and S. Tarantola, 2010,
\emph{Variance based sensitivity analysis of model output. Design and estimator for the
total sensitivity index}, Computer Physics Communications 181, 259--270.
}

\author{
Bertrand Iooss
}

\seealso{
  \code{\link{sobol}, \link{sobol2002}, \link{soboljansen}, \link{sobolEff}, \link{sobolmara}}
}

\examples{
# Test case : the non-monotonic Sobol g-function

# The method of sobol requires 2 samples
# There are 8 factors, all following the uniform distribution
# on [0,1]

library(boot)
n <- 1000
X1 <- data.frame(matrix(runif(8 * n), nrow = n))
X2 <- data.frame(matrix(runif(8 * n), nrow = n))

# sensitivity analysis

x <- sobol2007(model = sobol.fun, X1, X2, nboot = 100)
print(x)
plot(x)
}

\keyword{design}
