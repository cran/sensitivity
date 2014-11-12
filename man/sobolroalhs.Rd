\name{sobolroalhs}
\alias{sobolroalhs}
\alias{tell.sobolroalhs}
\alias{print.sobolroalhs}
\alias{plot.sobolroalhs}


\title{Sobol' Indices Estimation Using Replicated OA-based LHS}

\description{\code{sobolroalhs} implements the estimation of the Sobol' sensitivity indices introduced by Tissot & Prieur (2012) using two Orthogonal Array-based Latin Hypercubes. This function allows the estimation of all first-order indices or all closed second-order indices (containing the sum of the second-order effect between two inputs and the individual effects of each input) at a total cost of \eqn{2 \times N}{2*N}. For closed second-order indices \eqn{N=q^{2}}{N=q^2} where \eqn{q \geq d-1}{q >= d-1} is a prime number denoting the number of levels of the orthogonal array, and where \eqn{d}{d} is the number of factors.}

\usage{
sobolroalhs(model = NULL, factors, runs, order, conf=0.95, tail=TRUE, na.rm=FALSE, \dots)
\method{tell}{sobolroalhs}(x, y = NULL, \dots)
\method{print}{sobolroalhs}(x, \dots)
\method{plot}{sobolroalhs}(x, ylim = c(0,1), type="standard", \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{factors}{an integer specifying the number of factors, or a vector of
    character strings giving their names.}
  \item{runs}{an integer specifying the number \eqn{N}{N} of model runs.}
  \item{order}{an integer specifying the order of the indices (1 or 2).}
  \item{conf}{the confidence level for confidence intervals.}
  \item{tail}{a boolean specifying the method used to choose the number of levels of the orthogonal array (see first Warning messages).}
  \item{na.rm}{a boolean specifying if the response of the model contains \code{NA} values.}
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
  \item{OA}{the orthogonal array constructed (\code{NULL} if \code{order}=1).}
  \item{levels}{the number of levels of the orthogonal array constructed (\code{NULL} if \code{order}=1).}
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
    \item{"The number of model evaluations (runs) you entered is not the square of a prime number. It has been replaced by : "}{when \code{order}\eqn{=2}, the number of levels of the orthogonal array must be a prime number. If the number of runs specified is not a square of a prime number then this warning message indicates that the number of runs was replaced depending on the value of \code{tail}. If \code{tail=TRUE} (resp. \code{tail=FALSE}) the new number of runs is equals to the square of the prime number preceding (resp. following) the square root of \code{runs}.}
    
    \item{"The number of model evaluations (runs) you entered is not satisfying the constraint \eqn{n \geq (d-1)^2}{n >= (d-1)^2}. It has been replaced by : "}{when \code{order}\eqn{=2}, the number of runs must satisfied the constraint \eqn{N \geq (d-1)^{2}}{N \ge (d-1)^2} where \eqn{d}{d} is the number of factors. This warning message indicates that the number of runs was replaced by the square of the prime number following (or equals to) \eqn{d-1}.}}

}

\references{
Tissot, J. Y. and Prieur, C. (2012), \emph{Estimating Sobol's indices combining Monte Carlo integration and Latin hypercube sampling} \url{http://hal.archives-ouvertes.fr/hal-00743964}.

A.S. Hedayat, N.J.A. Sloane, John Stufken (1999), \emph{Orthogonal Arrays: Theory and Applications}.
}


\examples{
library(numbers)

# Test case : the non-monotonic Sobol g-function

# The method of sobol requires 2 samples
# (there are 8 factors, all following the uniform distribution on [0,1])

# first-order sensitivity indices
x <- sobolroalhs(model = sobol.fun, factors = 8, runs = 1000, order = 1)
print(x)
plot(x)

# closed second-order sensitivity indices
x <- sobolroalhs(model = sobol.fun, factors = 8, runs = 1000, order = 2)
print(x)
plot(x)

# Test case : the Ishigami function

# New function because sobolroalhs() works with U[0,1] inputs
ishigami1.fun=function(x) ishigami.fun(x*2*pi-pi)

# first-order sensitivity indices
x <- sobolroalhs(model = ishigami1.fun, factors = 3, runs = 100000, order = 1)
print(x)
plot(x)

# closed second-order sensitivity indices
x <- sobolroalhs(model = ishigami1.fun, factors = 3, runs = 100000, order = 2)
print(x)
plot(x)

# dealing with NA values
x <- sobolroalhs(model = NULL, factors = 3, runs = 100000, order =1,na.rm=TRUE)
y <- ishigami1.fun(x$X)
# we randomly insert NA values in y
pos <- sample(length(y),100)
y[pos] <- NA
tell(x,y)
print(x)
plot(x)
}
\keyword{design}

