\name{sobolSalt}
\alias{sobolSalt}
\alias{tell.sobolSalt}
\alias{print.sobolSalt}
\alias{plot.sobolSalt}
\alias{ggplot.sobolSalt}

\title{Monte Carlo Estimation of Sobol' Indices based on Saltelli schemes}

\description{
  \code{sobolSalt} implements the Monte Carlo estimation of
  the Sobol' indices for either both first-order and total effect indices at the same
  time (alltogether \eqn{2p}{2p} indices) at a total cost of \eqn{n\times(p+2)}{n*(p + 2)} model evaluations; or first-order, second-order and total indices at the same time (alltogether \eqn{2p+ p\times(p-1)/2}{2p+ p*(p-1)/2} indices) at a total cost of \eqn{n\times(2\times p+2)}{n*(2*p + 2)} model evaluations.
}

\usage{
sobolSalt(model = NULL, X1, X2, scheme="A", nboot = 0, conf = 0.95, \dots)
\method{tell}{sobolSalt}(x, y = NULL, \dots)
\method{print}{sobolSalt}(x, \dots)
\method{plot}{sobolSalt}(x, ylim = c(0, 1), choice, \dots)
\method{ggplot}{sobolSalt}(x, ylim = c(0, 1), choice, \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{X1}{the first random sample (containing \code{n} points).}
  \item{X2}{the second random sample (containing \code{n} points).}
  \item{scheme}{a letter \code{"A"} or \code{"B"} indicating which scheme to use (see "Details")}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level for bootstrap confidence intervals.}
  \item{x}{a list of class \code{"sobol"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{choice}{an integer specifying which indices to plot: \code{1} for first-order and total effect indices, \code{2} for second-order indices.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called}
}

\value{
  \code{sobolSalt} returns a list of class \code{"sobolSalt"}, containing all
  the input arguments detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments.}
  \item{y}{the response used.}
  \item{V}{the model variance.}
  \item{S}{the estimations of the Sobol' first-order indices.}
  \item{S2}{the estimations of the Sobol' second-order indices (only for scheme \code{"B"}).}
  \item{T}{the estimations of the Sobol' total sensitivity indices.}
}

\details{
The estimators used are the one implemented in \code{"sobolEff"}.

\code{scheme} specifies which Saltelli's scheme is to be used: \code{"A"} to estimate both first-order and total effect indices, \code{"B"} to estimate first-order, second-order and total effect indices.
}

\references{
A. Janon, T. Klein, A. Lagnoux, M. Nodet, C. Prieur (2014), \emph{Asymptotic normality and efficiency of two Sobol index estimators}, ESAIM: Probability and Statistics, 18:342-364.

A. Saltelli, 2002, \emph{Making best use of model evaluations to compute
sensitivity indices}, Computer Physics Communication, 145:580-297.
}

\author{
Laurent Gilquin
}

\seealso{
  \code{\link{sobol}, \link{sobol2007}, \link{soboljansen}, \link{sobolmartinez}, \link{sobolEff}}
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

x <- sobolSalt(model = sobol.fun, X1, X2, scheme="A", nboot = 100)
print(x)
plot(x, choice=1)

library(ggplot2)
ggplot(x, choice=1)
}

\keyword{design}
