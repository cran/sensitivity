\name{sobol2002}
\alias{sobol2002}
\alias{tell.sobol2002}
\alias{print.sobol2002}
\alias{plot.sobol2002}
\alias{plotMultOut.sobol2002}
\alias{ggplot.sobol2002}

\title{Monte Carlo Estimation of Sobol' Indices (scheme by Saltelli 2002)}

\description{
  \code{sobol2002} implements the Monte Carlo estimation of
  the Sobol' indices for both first-order and total indices at the same
  time (alltogether \eqn{2p}{2p} indices), at a total cost of \eqn{(p+2)
    \times n}{(p + 2) * n} model evaluations. These are called the Saltelli estimators.
}

\usage{
sobol2002(model = NULL, X1, X2, nboot = 0, conf = 0.95, \dots)
\method{tell}{sobol2002}(x, y = NULL, return.var = NULL, \dots)
\method{print}{sobol2002}(x, \dots)
\method{plot}{sobol2002}(x, ylim = c(0, 1), \dots)
  \method{plotMultOut}{sobol2002}(x, ylim = c(0, 1), \dots)
\method{ggplot}{sobol2002}(data, mapping = aes(), ylim = c(0, 1), \dots, environment
                 = parent.frame())
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{X1}{the first random sample.}
  \item{X2}{the second random sample.}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level for bootstrap confidence intervals.}
  \item{x}{a list of class \code{"sobol2002"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{data}{a list of class \code{"sobol2002"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{return.var}{a vector of character strings giving further
    internal variables names to store in  the output object \code{x}.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{mapping}{Default list of aesthetic mappings to use for plot. If not specified, 
    must be supplied in each layer added to the plot.}
  \item{environment}{[Deprecated] Used prior to tidy evaluation.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called}
}

\value{
  \code{sobol2002} returns a list of class \code{"sobol2002"}, containing all
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
BE CAREFUL! This estimator suffers from a conditioning problem when estimating 
the variances behind the indices computations. This can seriously affect the 
Sobol' indices estimates in case of largely non-centered output. To avoid this 
effect, you have to center the model output before applying \code{"sobol2002"}. 
Functions \code{"sobolEff"}, \code{"soboljansen"} and \code{"sobolmartinez"} 
do not suffer from this problem.
}

\references{
A. Saltelli, 2002, \emph{Making best use of model evaluations to compute
sensitivity indices}, Computer Physics Communication, 145, 580--297.
}

\author{
Gilles Pujol
}

\seealso{
  \code{\link{sobol}, \link{sobolSalt}, \link{sobol2007}, \link{soboljansen}, \link{sobolmartinez}, \link{sobolEff}, \link{sobolmara}, \link{sobolGP}, \link{sobolMultOut}}
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

x <- sobol2002(model = sobol.fun, X1, X2, nboot = 100)
print(x)
plot(x)

library(ggplot2)
ggplot(x)
}

\keyword{design}
