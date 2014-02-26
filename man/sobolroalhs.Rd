\name{sobolroalhs}
\alias{sobolroalhs}
\alias{tell.sobolroalhs}
\alias{print.sobolroalhs}
\alias{plot.sobolroalhs}


\title{Sobol' Indices Estimation Using Replicated OA-based LHS}

\description{\code{sobolroalhs} implements the estimation of the Sobol' sensitivity indices introduced by Tissot & Prieur (2012) using two Orthogonal Array-based Latin Hypercubes. This function allows the estimation of all first-order indices at a total cost of \eqn{2 \times N}{2*N} or all second-order subset indices (containing the sum of the second-order effect between two inputs and the individual effects of each input) at a cost of \eqn{2 \times q^{2}}{2*q^2} where \eqn{q}{q} is a prime number.}

\usage{
sobolroalhs(model = NULL, factors, levels, order, choice="A", conf=0.95, \dots)
\method{tell}{sobolroalhs}(x, y = NULL, \dots)
\method{print}{sobolroalhs}(x, \dots)
\method{plot}{sobolroalhs}(x, ylim = c(0, 1), type="standard", \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{factors}{an integer specifying the number of factors, or a vector of
    character strings giving their names.}
  \item{levels}{an integer specifying the number of levels of the design.}
  \item{order}{an integer specifying the order of the indices (1 or 2).}
  \item{choice}{a character (\code{A} or \code{B}) specifying the method to generate the orthogonal array-based design (\code{A} for the Bose method and \code{B} for the Hadamard generalized matrix based method).}
  \item{conf}{the confidence level for confidence intervals.}
  \item{x}{a list of class \code{"sobolroalhs"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{ylim}{coordinate plotting limits for the indices values.}
  \item{type}{a character specifying the type of estimator to plot (\code{standard} for the basic estimator or \code{monod} for the Janon-Monod estimator.)}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
}

\value{
  \code{sobolroalhs} returns a list of class \code{"sobolroalhs"}, containing all
  the input arguments detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{matrix} containing the design of experiments.}
  \item{y}{a vector of model responses.}
  \item{V}{a data.frame containing the estimations of the variance (\code{Vs} for the standard variance and \code{Veff} for the Janon-Monod variance).}
  \item{S}{a data.frame containing the estimations of the Sobol' sensitivity indices (\code{S} for the standard estimator and \code{Seff} for the Janon-Monod estimator).}

}

\details{
  The method used by sobolroalhs only considers models whose inputs follow uniform distributions on [0,1]. The transformations of each input (between its initial distribution and a U[0,1] distribution) have therefore to be realized before the call to sobolroalhs().

Bootstrap confidence intervals are not available with this method ; the given confidence intervals come from asymptotical formula.
}

\section{Warning messages}{
  \describe{
    \item{"the number of levels recquired was not a prime number. the number was replaced by : "}{when order=2, the number of levels must be a prime number. This warning message indicates that the number of levels specified was note a prime number and was automatically replaced by the prime number following the square root of \code{factors}.}
    
    \item{"the number of levels recquired was not satisfying the constraint, the number was replaced by : "}{when order=2, the number of levels must satisfied the constrain \eqn{N \geq (d-1)^{2}}{N \ge (d-1)^2} where \eqn{d}{d} is the number of factors.}}

}

\references{
Tissot, J. Y. and Prieur, C. (2012), \emph{Estimating Sobol's indices combining Monte Carlo integration and Latin hypercube sampling}.

A.S. Hedayat, N.J.A. Sloane, John Stufken (1999), \emph{Orthogonal Arrays: Theory and Applications}.
}


\seealso{
  \code{\link{sobol}, \link{sobolmara}}
}

\examples{
# Test case : the non-monotonic Sobol g-function
library(numbers)

# The method of sobol requires 2 samples
# (there are 8 factors, all following the uniform distribution on [0,1])

# first-order sensitivity indices
x <- sobolroalhs(model = sobol.fun, factors = 8, levels = 1000, order = 1)
print(x)
plot(x)

# global second-order sensitivity indices
x <- sobolroalhs(model = sobol.fun, factors = 8, levels = 1000, order = 2)
print(x)
plot(x)

# Test case : the Ishigami function

# New function because sobolroalhs() works with U[0,1] inputs
ishigami1.fun=function(x) ishigami.fun(x*2*pi-pi)

x <- sobolroalhs(model = ishigami1.fun, factors = 3, levels = 100000, order = 1)
print(x)
plot(x)

x <- sobolroalhs(model = ishigami1.fun, factors = 3, levels = 100000, order = 2)
print(x)
plot(x)
}

\keyword{design}

