\name{sobol}
\alias{sobol}
\alias{compute.sobol}
\alias{print.sobol}
\alias{plot.sobol}

\title{Sobol Non-linear Sensitivity Analysis}

\description{
  \code{sobol} is the implementation of the standard Monte Carlo
  estimation of the Sobol indices. This function generates the design of
  experiments based on two independent samples and computes the
  requested indices. These indices are ALL the Sobol anova-like
  decomposition indices, from the first order to the order given by the
  user.
}

\usage{
sobol(model = NULL, x1, x2, order = 1,
      nboot = 0, conf = 0.95, \dots)
\method{compute}{sobol}(sa, y = NULL)
}

\arguments{
  \item{model}{the model.}
  \item{x1}{the first random sample for factors.}
  \item{x2}{the second random sample for factors.}
  \item{order}{the maximum order of indices to compute.}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level for bootstrap confidence intervals.}
  \item{sa}{the sensitivity analysis object.}
  \item{y}{the response.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
}

\details{
  \code{model} is a function or a predictor (a class with a
  \code{predict} method) computing the response \code{y} based on the
  sample given by \code{x}. If no model is specified, the indices will be
  computed when one gives the response.

  The Monte Carlo estimation requires two independent random samples
  \code{x1} and \code{x2}. They must have the same dimensions.

  With this method, the estimation of the all the indices from the order 1 to
  the order \eqn{k}{k} (given by the \code{order} paremeter) require
  \deqn{n\sum_{j=0}^k \frac{p!}{j!(p-j)!}}{n * sum_{from j = 0 to k} p!
    / (j! * (p - j)!)} calls to the model, with \eqn{n}{n} the length of
  the samples \code{x1} and \code{x2}, and \eqn{p}{p} the number of factors.
}

\value{
  \code{sobol} returns an object of class \code{"sobol"}.
  An object of class \code{"sobol"} is a list containing the following
  components:

  \item{model}{the model.}
  \item{x1}{the first random sample for factors.}
  \item{x2}{the second random sample for factors.}
  \item{order}{the maximum order of indices to compute.}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level for bootstrap confidence intervals.}
  \item{x}{the factor sample.}
  \item{y}{the response used.}
  \item{S}{the estimations of the anova decomposition sensitivity indices.}
  \item{call}{the matched call.}
}

\references{
Saltelli, A., Chan, K. and Scott, E. M., 2000, \emph{Sensitivity
  analysis}. Wiley.

Sobol, I. M., 1993, Sensitivity analysis for non-linear mathematical
model. \emph{Math. Modelling Comput. Exp.}, \bold{1}, 407--414.
}

\seealso{
  \code{\link{sensitivity}}
  \code{\link{compute}}
}

\examples{
# Test case : the non-monotonic Sobol g-function

# The method of sobol requires 2 samples
# There are 8 factors, all following the uniform distribution
# on [0,1]

n <- 1000
x1 <- data.frame(matrix(nr = n, nc = 8))
x2 <- data.frame(matrix(nr = n, nc = 8))
for (i in 1:8){
  x1[, i] <- runif(n)
  x2[, i] <- runif(n)
}

# sensitivity analysis

sa <- sobol(model = sobol.fun, x1 = x1, x2 = x2, order = 2, nboot = 100)
print(sa)
#plot(sa)
}

\keyword{misc}
