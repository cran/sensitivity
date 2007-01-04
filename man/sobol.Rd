\name{sobol}
\alias{sobol}
\alias{sobol.sobol93}
\alias{tell.sobol.sobol93}
\alias{print.sobol.sobol93}
\alias{plot.sobol.sobol93}
\alias{sobol.saltelli02}
\alias{tell.sobol.saltelli02}
\alias{print.sobol.saltelli02}
\alias{plot.sobol.saltelli02}

\title{Sobol Non-linear Sensitivity Analysis}

\description{
  \code{sobol} is the implementation of the Monte Carlo estimation of
  the Sobol indices.
}

\usage{
sobol(method = "sobol93", model = NULL, x1, x2,
      max.order = 1, nboot = 0, conf = 0.95, \dots)
}

\arguments{
  \item{method}{the method: "sobol93" or "saltelli02"}
  \item{model}{the model}
  \item{x1}{the first random sample}
  \item{x2}{the second random sample}
  \item{max.order}{the maximum order of indices to compute (method "sobol93")}
  \item{nboot}{the number of bootstrap replicates}
  \item{conf}{the confidence level for bootstrap confidence intervals}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called}
}

\details{
  Two methods. The method "sobol93" computes all the Sobol indices
  (coming from the HDMR-ANOVA decomposition) from order 1 to order
  given by the argument \code{max.order}. The method "saltelli02"
  computes both first order and total indices at a reduced computational
  cost.
  
  \code{model} is a function or a predictor (a class with a
  \code{predict} method) computing the response \code{y} based on the
  sample given by \code{x}. If no model is specified, the indices will be
  computed when one gives the response.

  The Monte Carlo estimation requires two independent random samples
  \code{x1} and \code{x2}. They must have the same dimensions.
}

\value{
  \code{sobol} returns an object of class \code{"sobol"}.
  An object of class \code{"sobol"} is a list containing the following
  components:

  \item{x}{the factor sample}
  \item{y}{the response used}
  \item{S}{the estimations of the Sobol sensitivity indices}
  \item{St}{the estimations of the total indices (method "saltelli02")}
}

\section{Computational cost}{
  For the method "sobol93", the number of model evaluations is \eqn{n
    \times (N+1)}{n * (N + 1)} where \eqn{n}{n} is the size of the
  samples \code{x1} and \code{x2}, and \eqn{N}{N} is the number of indices
  to estimate.

  For the method "saltelli02", the number of model evaluations is \eqn{n
  \times (p+2)}{n * (p + 2)} where \eqn{p}{p} is the number of
  factors (for the estimation of \eqn{2p}{2p} indices).
}

\references{
Saltelli, A., Chan, K. and Scott, E. M., 2000, \emph{Sensitivity
  analysis}. Wiley.

Sobol, I. M., 1993, Sensitivity analysis for non-linear mathematical
model. \emph{Math. Modelling Comput. Exp.}, \bold{1}, 407--414.

Saltelli, A., 2002, Making best use of model evaluations to compute
sensitivity indices. \emph{Computer Physics Communication}, \bold{145},
580--297.
}

\examples{
# Test case : the non-monotonic Sobol g-function

# The method of sobol requires 2 samples
# There are 8 factors, all following the uniform distribution
# on [0,1]

n <- 1000
x <- data.frame(matrix(nr = 2 * n, nc = 8))
for (i in 1:8)
  x[, i] <- runif(2 * n)

# sensitivity analysis

sa <- sobol(model = sobol.fun, x1 = x[1:n,], x2 = x[(n+1):(2*n),], max.order = 2, nboot = 100)
print(sa)
#plot(sa)
}

\keyword{misc}
