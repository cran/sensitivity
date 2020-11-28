\name{sobolrank}
\alias{sobolrank}
\alias{tell.sobolrank}
\alias{print.sobolrank}
\alias{plot.sobolrank}
\alias{ggplot.sobolrank}

\title{First-order sensitivity indices estimation via ranking}

\description{ \code{sobolrank} implements the estimation of all first-order indices using only N model evaluations 
via ranking following Gamboa et al. (2020) and inspired by Chatterjee (2019).
}

\usage{
sobolrank(model = NULL, X, nboot = 0, conf = 0.95, nsample = round(0.8*nrow(X)), \dots)
\method{tell}{sobolrank}(x, y = NULL, \dots)
\method{print}{sobolrank}(x, \dots)
\method{plot}{sobolrank}(x, ylim = c(0, 1), \dots)
\method{ggplot}{sobolrank}(x, ylim = c(0, 1), \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{X}{a random sample of the inputs.}
  \item{nboot}{the number of bootstrap replicates, see details.}
  \item{conf}{the confidence level for confidence intervals, see details.}
  \item{nsample}{the size of the bootstrap sample, see details.}
  \item{x}{a list of class \code{"sobolrank"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
}

\value{
  \code{sobolrank} returns a list of class \code{"sobolrank"}, containing all
  the input arguments detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments.}
  \item{y}{a vector of model responses.}
  \item{S}{the estimations of the Sobol' sensitivity indices.}
}
\details{
		 The estimator used by sobolrank is defined in Gamboa et al. (2020).
     It is based on ranking the inputs as was first proposed by Chatterjee (2019) for a Cramer-Von Mises based estimator.
     All first-order indices can be estimated with a single sample of size N.
     Since boostrap creates ties which are not accounted for in the algorithm, confidence intervals are obtained by 
     sampling without replacement with a sample size \code{nsample}.
}
\references{
		 Gamboa, F., Gremaud, P., Klein, T., & Lagnoux, A. (2020). Global Sensitivity Analysis: 
		 a new generation of mighty estimators based on rank statistics. arXiv preprint arXiv:2003.01772.
		 
		 Chatterjee, S. (2019). A new coefficient of correlation. arXiv preprint arXiv:1909.10140.
}
\author{
Sebastien Da Veiga
}
\seealso{
  \code{\link{sobol}, \link{sobol2002}, \link{sobolSalt}, \link{sobol2007}, \link{soboljansen}, \link{sobolmartinez},
  \link{sobolSmthSpl}, \link{sobolEff}, \link{sobolshap_knn}}
}
\examples{
# Test case : the non-monotonic Sobol g-function
# Example with a call to a numerical model
library(boot)
n <- 1000
X <- data.frame(matrix(runif(8 * n), nrow = n))
x <- sobolrank(model = sobol.fun, X = X, nboot = 100)
print(x)
library(ggplot2)
ggplot(x)
# Test case : the Ishigami function
# Example with given data
n <- 500
X <- data.frame(matrix(-pi+2*pi*runif(3 * n), nrow = n))
Y <- ishigami.fun(X)
x <- sobolrank(model = NULL, X)
tell(x,Y)
print(x)
ggplot(x)
}