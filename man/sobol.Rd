\name{sobol}
\alias{sobol}
\alias{tell.sobol}
\alias{print.sobol}
\alias{plot.sobol}

\title{Monte Carlo Estimation of Sobol' Indices}

\description{ \code{sobol} implements the Monte Carlo estimation of
   the Sobol' sensitivity indices (standard estimator). This method allows the estimation of
   the indices of the variance decomposition, sometimes referred to as
   functional ANOVA decomposition, up to a given order, at a total cost
   of \eqn{(N+1) \times n}{(N + 1) * n} where \eqn{N}{N} is the number
   of indices to estimate. This function allows also the estimation of
   the so-called subset (or group) indices, i.e. the first-order indices with respect to 
   single multidimensional inputs.}

\usage{
sobol(model = NULL, X1, X2, order = 1, nboot = 0, conf = 0.95, \dots)
\method{tell}{sobol}(x, y = NULL, return.var = NULL, \dots)
\method{print}{sobol}(x, \dots)
\method{plot}{sobol}(x, ylim = c(0, 1), \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{X1}{the first random sample.}
  \item{X2}{the second random sample.}
  \item{order}{either an integer, the maximum order in the ANOVA
    decomposition (all indices up to this order will be computed), or a
    list of numeric vectors, the multidimensional compounds
    of the wanted subset indices.}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level for bootstrap confidence intervals.}
  \item{x}{a list of class \code{"sobol"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{return.var}{a vector of character strings giving further
    internal variables names to store in  the output object \code{x}.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
}

\value{
  \code{sobol} returns a list of class \code{"sobol"}, containing all
  the input arguments detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments.}
  \item{y}{a vector of model responses.}
  \item{V}{the estimations of Variances of the Conditional Expectations
    (VCE) with respect to one factor or one group of factors.}
  \item{D}{the estimations of the terms of the ANOVA decomposition (not
  for subset indices).}
  \item{S}{the estimations of the Sobol' sensitivity indices (not for
    subset indices).}

  Users can ask more ouput variables with the argument
  \code{return.var} (for example, bootstrap outputs \code{V.boot},
  \code{D.boot} and \code{S.boot}).
}

\references{
I. M. Sobol, 1993, \emph{Sensitivity analysis for non-linear mathematical
model}, Math. Modelling Comput. Exp., 1, 407--414.
}

\author{
Gilles Pujol
}

\seealso{
  \code{\link{sobol2002}, \link{sobolSalt}, \link{sobol2007}, \link{soboljansen}, \link{sobolmartinez},
  \link{sobolEff}, \link{sobolSmthSpl}, \link{sobolmara}, \link{sobolroalhs}, \link{fast99}, \link{sobolGP},
  \link{sobolMultOut}}
}

\examples{
# Test case : the non-monotonic Sobol g-function

# The method of sobol requires 2 samples
# (there are 8 factors, all following the uniform distribution on [0,1])
library(boot)
n <- 1000
X1 <- data.frame(matrix(runif(8 * n), nrow = n))
X2 <- data.frame(matrix(runif(8 * n), nrow = n))

# sensitivity analysis
x <- sobol(model = sobol.fun, X1 = X1, X2 = X2, order = 2, nboot = 100)
print(x)
#plot(x)
}

\keyword{design}
