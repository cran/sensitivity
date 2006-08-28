\name{fast}
\alias{fast}
\alias{compute.fast}
\alias{print.fast}
\alias{plot.fast}

\title{Fourier Amplitude Sensitivity Test}

\description{
  \code{fast} is the implementation of the Fourier Amplitude Sensitivity
  Test. This function computes the following non-linear indices of
  sensitivity: first-order indices and total indices.
}

\usage{
fast(model = NULL, factors, n, M = 4,
     G = "uniform", min = 0, max = 1, omega = NULL, \dots)
\method{compute}{fast}(sa, y = NULL)
}

\arguments{
  \item{model}{the model.}
  \item{factors}{the number of factors, or their names.}
  \item{n}{the length of the discretization of the s-space.}
  \item{M}{the number of harmonics to sum.}
  \item{G}{the space transformation.}
  \item{min}{the minimum values for the factors (uniform distribution).}
  \item{max}{the maximum values for the factors (uniform distribution).}
  \item{omega}{the set of frequencies.}
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

  The space transformation \code{G} is the function such that:
  \deqn{
    x_i = G_i(\mathrm{sin}(\omega_i s))
  }{
    x_i = G_i( sin( omega_i * s ) )
  }
  It must be a function of two parameters (\code{G <- function(i, x)
    ...}). If the string \code{"uniform"} is given, then the function is
  the best one for uniform factors on the range \eqn{[a_i, b_i]}{[a_i,
    b_i]}:
  \deqn{
    G_i(x) = a_i + (b_i-a_i) \times \left( \frac{1}{2} + \frac{1}{\pi}
    \mathrm{asin}(x) \right)
  }{
    G_i(x) = a_i + (b_i - a_i) * ( 1/2 + 1/pi * asin(x) )
  }
  where \eqn{a_i}{a_i} and \eqn{b_i}{b_i} are the boundaries given by
  the arguments \code{min} and \code{max}. \code{min} and \code{max} can
  be single values (the same for each factor) or vectors.

  If the set of frequencies \code{omega} is not given, the function use
  the set recommended by Saltelli (the so-called extended-FAST). The
  first frequency is the greater, corresponding to the index of
  interest, and the other correspond to the complementary set.
}

\value{
  \code{fast} returns an object of class \code{"fast"}.
  An object of class \code{"fast"} is a list containing the following
  components:

  \item{model}{the model.}
  \item{M}{the number of harmonics to sum.}
  \item{s}{the discretisation of the s-space.}
  \item{omega}{the set of frequencies.}
  \item{x}{the factor sample.}
  \item{y}{the response.}
  \item{S}{the estimations of the first-order indices.}
  \item{St}{the estimations of the total indices.}
  \item{call}{the matched call.}
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

\seealso{
  \code{\link{sensitivity}}
  \code{\link{compute}}
}

\examples{
# Test case : the non-monotonic Sobol g-function

sa <- fast(model = sobol.fun, factors = 8, n = 1000)
print(sa)
plot(sa)
}

\keyword{misc}
