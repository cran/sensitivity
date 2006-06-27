\name{morris}
\alias{morris}
\alias{compute.morris}
\alias{print.morris}
\alias{plot.morris}

\title{The Morris OAT Screening Method}

\description{
  \code{morris} is the implementation of the Morris OAT Screening
  method. This function generates the Morris design of experiments and
  computes the measures of sensitivity \eqn{\mu^*}{mu*} and
  \eqn{\sigma}{sigma}.
}

\usage{
morris(model = NULL, factors, levels, r, k.delta = "usual",
       min = 0, max = 1, scale = TRUE, nboot = 0, conf = 0.95, \dots)
\method{compute}{morris}(sa, y = NULL)
}

\arguments{
  \item{model}{the model.}
  \item{factors}{the number of factors, or their names.}
  \item{levels}{the number of levels of the design grid.}
  \item{r}{the number of repetitions of the design, i.e. the number of
    elementary effect computed per factor.}
  \item{k.delta}{the \sQuote{grid jump} coefficient.}
  \item{min}{the minimum values for the factors.}
  \item{max}{the maximum values for the factors.}
  \item{scale}{logical. If \code{TRUE}, the input and output data are
    scaled.}
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

  The number of levels is the same for each space coordinate. Then
  \code{levels} must be a single integer.

  \code{k.delta} is such that:
  \deqn{
    \Delta = \frac{\textrm{k.delta}}{k-1}
  }{
    Delta = \code{k.delta} / ( k - 1 )
  }
  where \eqn{k}{k} is the number of levels (\code{levels}). If
  \code{k.delta} is given as \code{"usual"} and \eqn{k}{k} is even,
  then \eqn{\Delta}{Delta} is the value recommended by Morris:
  \deqn{
    \Delta = \frac{k}{2(k-1)}
  }{
    Delta = k / ( 2 * ( k - 1 ) )
  }

  \code{min} and \code{max} are boundaries of the region of
  experimentation. They can be single values (the same for each
  factor) or vectors.
}

\value{
  \code{morris} returns an object of class \code{"morris"}.
  An object of class \code{"morris"} is a list containing the following
  components:

  \item{model}{the model.}
  \item{levels}{the number of levels of the design grid.}
  \item{r}{the number of repetitions of the design.}
  \item{delta}{the value of \eqn{\Delta}{Delta}.}
  \item{min}{the minimum values for the factors.}
  \item{max}{the maximum values for the factors.}
  \item{scale}{logical. If \code{TRUE}, the input and output data are
    scaled before computing the elementary effects.}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level for bootstrap confidence intervals.}
  \item{D}{the successive diagonal matrices composed of equiprobable +1 and -1.}
  \item{P}{the successive random permutation matrices.}
  \item{x}{the design of experiments (input sample).}
  \item{y}{the response.}
  \item{mu}{the estimations of the \eqn{\mu^*}{mu*} index.}
  \item{sigma}{the estimations of the \eqn{\sigma}{sigma} index.}
  \item{call}{the matched call.}
}

\references{
Saltelli, A., Chan, K. and Scott, E. M., 2000, \emph{Sensitivity
  analysis}. Wiley.

Morris, M. D., 1991, Factorial sampling plans for preliminary
computational experiments. \emph{Technometrics}, \bold{33}, 161--174.
}

\seealso{
  \code{\link{sensitivity}}
  \code{\link{compute}}
}

\examples{
# Test case : the non-monotonic function of Morris

sa <- morris(model = morris.fun, factors = 20, levels = 4, r = 4)
print(sa)
plot(sa)
}

\keyword{misc}
