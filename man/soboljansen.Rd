\name{soboljansen}
\alias{soboljansen}
\alias{tell.soboljansen}
\alias{print.soboljansen}
\alias{plot.soboljansen}

\title{Monte Carlo Estimation of Sobol' Indices (improved formulas of Jansen (1999) and Saltelli et al. (2010))}

\description{
  \code{soboljansen} implements the Monte Carlo estimation of
  the Sobol' indices for both first-order and total indices at the same
  time (alltogether \eqn{2p}{2p} indices), at a total cost of \eqn{(p+2)
    \times n}{(p + 2) * n} model evaluations. These are called the Jansen estimators.
}

\usage{
soboljansen(model = NULL, X1, X2, nboot = 0, conf = 0.95, \dots)
\method{tell}{soboljansen}(x, y = NULL, return.var = NULL, \dots)
\method{print}{soboljansen}(x, \dots)
\method{plot}{soboljansen}(x, ylim = c(0, 1), \dots)
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
  \code{soboljansen} returns a list of class \code{"soboljansen"}, containing all
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
This estimator is good for large first-order indices, and (large and small) total indices.
}

\references{
M.J.W. Jansen, 1999, \emph{Analysis of variance designs for model output},
Computer Physics Communication, 117, 35--43.

A. Saltelli, P. Annoni, I. Azzini, F. Campolongo, M. Ratto and S. Tarantola, 2010,
\emph{Variance based sensitivity analysis of model output. Design and estimator for the
total sensitivity index}, Computer Physics Communications 181, 259--270.
}

\author{
Bertrand Iooss
}

\seealso{
  \code{\link{sobol}, \link{sobol2002}, \link{sobol2007}, \link{sobolmartinez}, \link{sobolEff}, \link{sobolmara}}
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

x <- soboljansen(model = sobol.fun, X1, X2, nboot = 100)
print(x)
plot(x)
}

\keyword{design}
