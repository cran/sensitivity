\name{sobolroalhs}
\alias{sobolroalhs}
\alias{tell.sobolroalhs}
\alias{print.sobolroalhs}
\alias{plot.sobolroalhs}


\title{Sobol' Indices Estimation Using Replicated OA-based LHS}

\description{\code{sobolroalhs} implements the estimation of the Sobol' sensitivity indices introduced by Tissot & Prieur (2015) using two replicated designs (Latin hypercubes or orthogonal arrays). This function estimates either all first-order indices or all closed second-order indices at a total cost of \eqn{2 \times N}{2*N} model evaluations. For closed second-order indices \eqn{N=q^{2}}{N=q^2} where \eqn{q \geq d-1}{q >= d-1} is a prime number corresponding to the number of levels of the orthogonal array, and where \eqn{d}{d} indicates the number of factors.}

\usage{
sobolroalhs(model = NULL, factors, N, p=1, order, tail=TRUE, conf=0.95, nboot=0, \dots)
\method{tell}{sobolroalhs}(x, y = NULL, \dots)
\method{print}{sobolroalhs}(x, \dots)
\method{plot}{sobolroalhs}(x, ylim = c(0,1), \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{factors}{an integer giving the number of factors, or a vector of character strings giving their names.}
  \item{N}{an integer giving the size of each replicated design (for a total of \eqn{2 \times N}{2*N} model evaluations).}
  \item{p}{an integer giving the number of model outputs.}
  \item{order}{an integer giving the order of the indices (1 or 2).}
  \item{tail}{a boolean specifying the method used to choose the number of levels of the orthogonal array (see "Warning messages").}
  \item{conf}{the confidence level for confidence intervals.}
  \item{nboot}{the number of bootstrap replicates.}
  \item{x}{a list of class \code{"sobolroalhs"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
}

\value{
  \code{sobolroalhs} returns a list of class \code{"sobolroalhs"}, containing all
  the input arguments detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments (row concatenation of the two replicated designs).}
  \item{y}{the responses used.}
  \item{OA}{the orthogonal array constructed (\code{NULL} if \code{order}=1).}
  \item{V}{the estimations of Variances of the Conditional Expectations (VCE) with respect to each factor.}
  \item{S}{the estimations of the Sobol' indices.}

}

\details{

\code{sobolroalhs} automatically assigns a uniform distribution on [0,1] to each input. Transformations of distributions (between U[0,1] and the wanted distribution) have to be realized before the call to tell() (see "Examples").

Missing values (i.e \code{NA} values) in outputs are automatically handled by the function.

This function also supports multidimensional outputs (matrices in \code{y} or as output of \code{model}).
In this case, aggregated Sobol' indices are returned (\code{see sobolMultOut}).
}

\section{Warning messages}{
  \describe{
    \item{"The value entered for \code{N} is not the square of a prime number. It has been replaced by: "}{when \code{order}\eqn{=2}, the number of levels of the orthogonal array must be a prime number. If \code{N} is not a square of a prime number, then this warning message indicates that it was replaced depending on the value of \code{tail}. If \code{tail=TRUE} (resp. \code{tail=FALSE}) the new value of \code{N} is equal to the square of the prime number preceding (resp. following) the square root of \code{N}.}
    
    \item{"The value entered for \code{N} is not satisfying the constraint \eqn{N \geq (d-1)^2}{N >= (d-1)^2}. It has been replaced by: "}{when \code{order}\eqn{=2}, the following constraint must be satisfied \eqn{N \geq (d-1)^{2}}{N \ge (d-1)^2} where \eqn{d}{d} is the number of factors. This warning message indicates that \code{N} was replaced by the square of the prime number following (or equals to) \eqn{d-1}.}}

}

\references{
A.S. Hedayat, N.J.A. Sloane and J. Stufken, 1999, \emph{Orthogonal Arrays: Theory and Applications}, Springer Series in Statistics.

F. Gamboa, A. Janon, T. Klein and A. Lagnoux, 2014, \emph{Sensitivity indices for multivariate outputs}, Electronic Journal of Statistics, 8:575-603.

J.Y. Tissot and C. Prieur, 2015, \emph{Estimating Sobol's indices combining Monte Carlo integration and Latin hypercube sampling}, J. Statist. Comput. Simulation, 85:1358-1381.
}

\author{
Laurent Gilquin
}

\seealso{
  \code{\link{sobolmara}},
  \code{\link{sobolroauc}},
  \code{\link{sobolMultOut}}
}

\examples{
library(boot)
library(numbers)

####################
# Test case: the non-monotonic Sobol g-function

# The method of sobol requires 2 samples
# (there are 8 factors, all following the uniform distribution on [0,1])

# first-order sensitivity indices
x <- sobolroalhs(model = sobol.fun, factors = 8, N = 1000, order = 1, nboot=100)
print(x)
plot(x)

# closed second-order sensitivity indices
x <- sobolroalhs(model = sobol.fun, factors = 8, N = 1000, order = 2, nboot=100)
print(x)
plot(x)

####################
# Test case: dealing with non-uniform distributions

x <- sobolroalhs(model = NULL, factors = 3, N = 1000, order =1, nboot=0)

# X1 follows a log-normal distribution:
x$X[,1] <- qlnorm(x$X[,1])

# X2 follows a standard normal distribution:
x$X[,2] <- qnorm(x$X[,2])

# X3 follows a gamma distribution:
x$X[,3] <- qgamma(x$X[,3],shape=0.5)

# toy example
toy <- function(x){rowSums(x)}
y <- toy(x$X)
tell(x, y)
print(x)
plot(x)

####################
# Test case : multidimensional outputs

toy <- function(x){cbind(x[,1]+x[,2]+x[,1]*x[,2],2*x[,1]+3*x[,1]*x[,2]+x[,2])}
x <- sobolroalhs(model = toy, factors = 3, N = 1000, p=2, order =1, nboot=100)
print(x)
plot(x)
}
\keyword{design}

