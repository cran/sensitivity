\name{testmodels}
\alias{testmodels}
\alias{sobol.fun}
\alias{ishigami.fun}
\alias{morris.fun}

\title{Test Models For Sensitivity Analysis}

\description{
  These functions are standard testcase for sensitivity
  analysis benchmarks. There are: the g-function of Sobol, the function
  of Ishigami and the function of Morris.
}

\usage{
sobol.fun(x)
ishigami.fun(x)
morris.fun(x)
}

\arguments{
  \item{x}{the \code{matrix} or \code{data.frame} containing the input
    values.}
}

\value{
  All these functions return a \code{numeric} vector containig the values
  of the function. These functions are vectorized.
}

\references{
Saltelli, A., Chan, K. and Scott, E. M., 2000, \emph{Sensitivity
  analysis}. Wiley, 33--45.
}

\keyword{misc}
