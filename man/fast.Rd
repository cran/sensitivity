\name{fast}
\alias{fast}
\alias{fast.saltelli99}
\alias{tell.fast.saltelli99}
\alias{print.fast.saltelli99}
\alias{plot.fast.saltelli99}

\title{Fourier Amplitude Sensitivity Test}

\description{
  \code{fast} is the implementation of the Fourier Amplitude Sensitivity
  Test.
}

\usage{
fast(method = "saltelli99", model = NULL, factors, n, M = 4,
     omega = NULL, q = NULL, q.arg = NULL, \dots)
}

\arguments{
  \item{method}{the method: "saltelli99" only}
  \item{model}{the model}
  \item{factors}{the number of factors, or their names}
  \item{n}{the sample size}
  \item{M}{the interference parameter}
  \item{omega}{the set of frequencies}
  \item{q}{the names of the quantile functions for the factors distributions}
  \item{q.arg}{the quantile parameters}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called}
}

\details{
  The method "saltelli99" is the so-called extended fast method wich
  provides estimations of both first order and total indices at a
  low computational cost.
  
  \code{model} is a function or a predictor (a class with a
  \code{predict} method) computing the response \code{y} based on the
  sample given by \code{x}. If no model is specified, the indices will be
  computed when one gives the response.

  \code{factors} could either be a single number or a vector of
  character strings.
  
  \code{n} is the length of the discretization of the s-space (for
  computing Fourier coefficients) and \code{M} is the number of
  harmonics to sum (for computing partial variances).
  
  If the set of frequencies \code{omega} is not given, the function use
  the set recommended by the corresponding method. For the method
  "saltelli99", the first frequency is the greater, associated with the
  input variable to assess, and the other frequencies are associated
  with the complementary set.

  If \code{q} and \code{q.args} are not given, the factors will be
  considerd uniform on [0,1]. \code{q} is a list of character strings
  giving the names of the quantile functions (one for each factor), such
  as \code{qunif}, \code{qnorm}\dots \code{q} could  also be a single
  character string (the same for all). \code{q.arg} is a list of lists, each
  list being additional parameters for the corresponding quantile
  function. For example, the parameters of the quantile \code{qunif}
  could be \code{(min=1, max=2)} giving an uniform distribution on [1,2].
  If \code{q} is a single character string, then \code{q.arg} must
  be a single list.
}

\value{
  \code{fast} returns an object of class \code{"fast"}.
  An object of class \code{"fast"} is a list containing the following
  components:

  \item{x}{the factor sample}
  \item{y}{the response}
  \item{S}{the estimations of the first-order indices}
  \item{St}{the estimations of the total indices (method "saltelli99")}
}

\section{Computational cost}{
  For the method "saltelli99", the number of model evaluations is
  \eqn{p \times n}{p * n} where \eqn{p}{p} is the number of factors.
}

\references{
Saltelli, A., Tarantola, S. and Chan, K., 1999, A quantitative, model
independent method for global sensitivity analysis of model
output. \emph{Technometrics}, \bold{41}, 39--56.

Saltelli, A., Chan, K. and Scott, E. M., 2000, \emph{Sensitivity
  analysis}. Wiley.

Cukier, R. I., Levine, H. B. and Schuler, K. E., 1978, Nonlinear
sensitivity analysis of multiparameter model
systems. \emph{J. Comput. Phys.}, \bold{26}, 1--42.
}

\examples{
# Test case : the non-monotonic Ishigami function

sa <- fast(model = ishigami.fun, factors = 3, n = 1000,
           q = "qunif", q.arg = list(min = -pi, max = pi))
print(sa)
plot(sa)
}

\keyword{misc}
