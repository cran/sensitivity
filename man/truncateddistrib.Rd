\name{truncateddistrib}
\alias{truncateddistrib}
\alias{dnorm.trunc}
\alias{pnorm.trunc}
\alias{qnorm.trunc}
\alias{rnorm.trunc}
\alias{dgumbel.trunc}
\alias{pgumbel.trunc}
\alias{qgumbel.trunc}
\alias{rgumbel.trunc}

\title{Truncated distributions}

\description{
  \code{dnorm.trunc}, \code{pnorm.trunc}, \code{qnorm.trunc} and 
  \code{rnorm.trunc} are functions for the Truncated Normal Distribution.
  \code{dgumbel.trunc}, \code{pgumbel.trunc}, \code{qgumbel.trunc} and 
  \code{rgumbel.trunc} are functions for the Truncated Gumbel Distribution.
}

\usage{
dnorm.trunc(x, mean = 0, sd = 1, min = -1e6, max = 1e6)
pnorm.trunc(q, mean = 0, sd = 1, min = -1e6, max = 1e6)
qnorm.trunc(p, mean = 0, sd = 1, min = -1e6, max = 1e6)
rnorm.trunc(n, mean = 0, sd = 1, min = -1e6, max = 1e6)
dgumbel.trunc(x, loc = 0, scale = 1, min = -1e6, max = 1e6)
pgumbel.trunc(q, loc = 0, scale = 1, min = -1e6, max = 1e6)
qgumbel.trunc(p, loc = 0, scale = 1, min = -1e6, max = 1e6)
rgumbel.trunc(n, loc = 0, scale = 1, min = -1e6, max = 1e6)
}

\arguments{
  \item{x, q}{vector of quantiles}
  \item{p}{vector of probabilities}
  \item{n}{number of observations}
  \item{mean, sd}{means and standard deviation parameters}
  \item{loc, scale}{location and scale parameters}
  \item{min}{vector of minimal bound values}
  \item{max}{vector of maximal bound values}
}

\details{
  See \code{dnorm} for details on the Normal distribution.
  The Gumbel distribution comes from the evd package.
  See \code{dgumbel} for details on the Gumbel distribution.
}

\value{
  \code{dnorm.trunc} and \code{dgumbel.trunc} give the density, \code{pnorm} and \code{pgumbel.trunc} give the distribution function, \code{qnorm} and \code{qgumbel.trunc} give the quantile function, \code{rnorm} and \code{rgumbel.trunc} generate random deviates.
}

\author{
Gilles Pujol and Bertrand Iooss
}

\keyword{misc}
