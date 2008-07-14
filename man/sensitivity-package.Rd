\name{sensitivity-package}
\alias{sensitivity-package}
\alias{sensitivity}
\docType{package}

\title{Sensitivity Analysis}

\description{Methods and functions for global sensitivity analysis.}

\details{
  The \pkg{sensitivity} package implements some global sensitivity analysis
  methods:
  \itemize{
    \item Linear regression coefficients: SRC and SRRC
    (\code{\link{src}}), PCC and PRCC (\code{\link{pcc}}).
    \item Morris's "OAT" elementary effects screening method (\code{\link{morris}}).
    \item Bettonvil's sequential bifurcations (\code{\link{sb}}).
    \item Monte Carlo estimation of Sobol' indices: Sobol's scheme
  (1993) to compute the indices given by the variance decomposition up
  to a specified order (\code{\link{sobol}}), and Saltelli's
    scheme (2002) to compute first order and total indices
    with a reduced cost (\code{\link{sobol2002}}).
    \item Estimation of the Sobol' first order and total indices with
    Saltelli's so-called "extended-FAST" method (\code{\link{fast99}}).
  }
  
  Moreover, some utilities are provided: standard test-cases
  (\code{\link{testmodels}}) and template file generation
  (\code{\link{template.replace}}).
}

\section{Model managing}{
  The \pkg{sensitivity} package works either on \R models than on
  external models (such as executables).

  \R models must be functions or objects that have a \code{predict}
  method, such as \code{lm} objects.  Models defined as functions will
  be called once with an expression of the form \code{y <- f(X)} where
  \code{X} is the design of experiments, i.e. a \code{data.frame} with
  \code{p} columns (the input factors) and \code{n} lines (each, an
  experiment), and \code{y} is the vector of length \code{n} of the
  model responses (we say that such functions are vectorized).

  If the model is external to \R, for instance a computational code, it
  must be analyzed with the decoupled approach, see
  \code{\link{decoupling}}. This approach can also be used on \R models
  that doesn't fit the specifications.
}

\references{
A. Saltelli, K. Chan and E. M. Scott eds, 2000, \emph{Sensitivity
  Analysis}, Wiley.
}

\keyword{package}

