\name{fast99}
\alias{fast99}
\alias{tell.fast99}
\alias{print.fast99}
\alias{plot.fast99}

\title{Extended Fourier Amplitude Sensitivity Test}

\description{
  \code{fast99} implements the so-called "extended-FAST" method
  (Saltelli et al. 1999). This method allows the estimation of first
  order and total Sobol' indices for all the factors (alltogether
  \eqn{2p}{2p} indices, where \eqn{p}{p} is the number of factors) at a
  total cost of \eqn{n \times p}{n * p} simulations.
}

\usage{
fast99(model = NULL, factors, n, M = 4, omega = NULL,
       q = NULL, q.arg = NULL, \dots)
\method{tell}{fast99}(x, y = NULL, \dots)
\method{print}{fast99}(x, \dots)
\method{plot}{fast99}(x, ylim = c(0, 1), \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{factors}{an integer giving the number of factors, or a vector of
    character strings giving their names.}
  \item{n}{an integer giving the sample size, i.e. the length of the
  discretization of the s-space (see Cukier et al.).}
  \item{M}{an integer specifying the interference parameter, i.e. the
  number of harmonics to sum in the Fourier series decomposition (see
  Cukier et al.).}
  \item{omega}{a vector giving the set of frequencies, one frequency for
  each factor (see details below).}
  \item{q}{a vector of quantile functions names corresponding to
  wanted factors distributions (see details below).}
  \item{q.arg}{a list of quantile functions parameters (see details below).}
  \item{x}{a list of class \code{"fast99"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
}

\details{
  
  If not given, the set of frequencies \code{omega} is taken from
  Saltelli et al.  The first frequency of the vector \code{omega} is
  assigned to each factor \eqn{X_i}{X_i} in turn (corresponding to the
  estimation of Sobol' indices \eqn{S_i}{S_i} and \eqn{S_{T_i}}{ST_i}),
  other frequencies being assigned to the remaining factors.

  If the arguments \code{q} and \code{q.args} are not given, the factors
  are taken uniformly distributed on \eqn{[0,1]}{[0,1]}. The
  argument \code{q} must be list of character strings, giving the names
  of the quantile functions (one for each factor), such as \code{qunif},
  \code{qnorm}\dots It can also be a single character string, meaning
  same distribution for all. The argument \code{q.arg} must be a list of
  lists, each one being additional parameters for the corresponding
  quantile function. For example, the parameters of the quantile
  function \code{qunif} could be \code{list(min=1, max=2)}, giving an
  uniform distribution on \eqn{[1,2]}{[1,2]}.  If \code{q} is a single
  character string, then \code{q.arg} must be a single list (rather than
  a list of one list).
  
}

\value{
  \code{fast99} returns a list of class \code{"fast99"}, containing all
  the input arguments detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the factors sample values.}
  \item{y}{a vector of model responses.}
  \item{V}{the estimation of variance.}
  \item{D1}{the estimations of Variances of the Conditional Expectations
  (VCE) with respect to each factor.}
  \item{Dt}{the estimations of VCE with respect to each factor
  complementary set of factors ("all but \eqn{X_i}{Xi}").}
}


\references{
A. Saltelli, S. Tarantola and K. Chan, 1999, \emph{A quantitative, model
independent method for global sensitivity analysis of model
output}, Technometrics, 41, 39--56.

R. I. Cukier, H. B. Levine and K. E. Schuler, 1978, \emph{Nonlinear
sensitivity analysis of multiparameter model
systems}. J. Comput. Phys., 26, 1--42.
}

\seealso{
  \code{\link{decoupling}}
}

\examples{
# Test case : the non-monotonic Ishigami function
x <- fast99(model = ishigami.fun, factors = 3, n = 1000,
            q = "qunif", q.arg = list(min = -pi, max = pi))
print(x)
plot(x)
}

\keyword{design}
