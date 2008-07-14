\name{testmodels}
\alias{testmodels}
\alias{sobol.fun}
\alias{ishigami.fun}
\alias{morris.fun}

\title{Test Models for Sensitivity Analysis}

\description{
  These functions are standard testcase for sensitivity
  analysis benchmarks: The g-function of Sobol', the function
  of Ishigami and the function of Morris (see Saltelli et al. 2000,
  section 2.9)
}

\usage{
sobol.fun(X)
ishigami.fun(X)
morris.fun(X)
}

\arguments{
  \item{X}{a matrix (or \code{data.frame}) containing the input
    sample.}
}

\value{
  A vector of function responses.
}

\references{
A. Saltelli, K. Chan and E. M. Scott eds, 2000, \emph{Sensitivity
  Analysis}, Wiley.
}

\keyword{misc}
