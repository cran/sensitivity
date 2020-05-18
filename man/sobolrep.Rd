\name{sobolrep}
\alias{sobolrep}
\alias{tell.sobolrep}
\alias{print.sobolrep}
\alias{plot.sobolrep}


\title{Sobol' indices estimation based on replicated orthogonal arrays}

\description{\code{sobolrep} generalizes the estimation of the Sobol' sensitivity indices introduced by Tissot & Prieur (2015) using two replicated orthogonal arrays. This function estimates either 
\itemize{
\item all first-order and second-order indices at a total cost of \eqn{2 \times N}{2*N} model evaluations,
\item or all first-order, second-order and total-effect indices at a total cost of \eqn{N \times (d+2)}{N*(d+2)} model evaluations,}
where \eqn{N=q^{2}}{N=q^2} and \eqn{q \geq d-1}{q >= d-1} is a prime number corresponding to the number of levels of each orthogonal array.}

\usage{
sobolrep(model = NULL, factors, N, tail=TRUE, 
			conf=0.95, nboot=0, nbrep=1, total=FALSE, \dots)
\method{tell}{sobolrep}(x, y = NULL, \dots)
\method{print}{sobolrep}(x, \dots)
\method{plot}{sobolrep}(x, ylim = c(0,1), choice, \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{factors}{an integer giving the number of factors, or a vector of character strings giving their names.}
  \item{N}{an integer giving the size of each replicated design (for a total of \eqn{2 \times N}{2*N} model evaluations).}
    \item{tail}{a boolean specifying the method used to choose the number of levels of the orthogonal array (see "Warning messages").}
  \item{conf}{the confidence level for confidence intervals.}
  \item{nboot}{the number of bootstrap replicates.}
  \item{nbrep}{the number of times the estimation procedure is repeated (see "Details").}
  \item{total}{a boolean specifying whether or not total effect indices are estimated.}
  \item{x}{a list of class \code{"sobolrep"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{the model response.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{choice}{an integer specifying which indices to plot: \code{1} for first-order indices, \code{2} for second-order indices, \code{3} for total-effect indices.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
}

\value{
  \code{sobolrep} returns a list of class \code{"sobolrep"}, containing all
  the input arguments detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments (row concatenation of the two replicated designs).}
  \item{y}{the response used.}
  \item{RP}{the matrix of permutations.}
  \item{V}{the model variance.}
  \item{S}{a data.frame containing estimations of the first-order Sobol' indices plus confidence intervals if specified.}
  \item{S2}{a data.frame containing estimations of the second-order Sobol' indices plus confidence intervals if specified.}
  \item{T}{a data.frame containing estimations of the total-effect indices plus confidence intervals if specified.}

}

\details{

\code{sobolrep} automatically assigns a uniform distribution on [0,1] to each input. Transformations of distributions (between U[0,1] and the wanted distribution) have to be performed before the call to tell() (see "Examples").

\code{nbrep} specifies the number of times the estimation procedure is repeated. Each repetition makes use of the orthogonal array structure to obtain a new set of Sobol' indices. It is important to note that no additional model evaluations are performed (the cost of the procedure remains the same).
}

\section{Warning messages}{
  \describe{
    \item{"The value entered for \code{N} is not the square of a prime number. It has been replaced by: "}{the number of levels \code{q} of each orthogonal array must be a prime number. If \code{N} is not a square of a prime number, then this warning message indicates that it was replaced depending on the value of \code{tail}. If \code{tail=TRUE} (resp. \code{tail=FALSE}) the new value of \code{N} is equal to the square of the prime number preceding (resp. following) the square root of \code{N}.}
    
    \item{"The value entered for \code{N} is not satisfying the constraint \eqn{N \geq (d-1)^2}{N >= (d-1)^2}. It has been replaced by: "}{the following constraint must be satisfied \eqn{N \geq (d-1)^{2}}{N \ge (d-1)^2} where \eqn{d}{d} is the number of factors. This warning message indicates that \code{N} was replaced by the square of the prime number following (or equals to) \eqn{d-1}.}}

}

\references{
A.S. Hedayat, N.J.A. Sloane and J. Stufken, 1999, \emph{Orthogonal Arrays: Theory and Applications}, Springer Series in Statistics.

J-Y. Tissot and C. Prieur, 2015, \emph{A randomized orthogonal orray-based procedure for the estimation of first- and second-order Sobol' indices}, J. Statist. Comput. Simulation, 85:1358-1381.
}


\examples{
# Test case: the non-monotonic Sobol g-function

# The method of sobol requires 2 samples
# (there are 8 factors, all following the uniform distribution on [0,1])

x <- sobolrep(model = sobol.fun, factors = 8, N = 1000, nboot=100, nbrep=1, total=FALSE)
print(x)
plot(x,choice=1)
plot(x,choice=2)

# Test case: dealing with non-uniform distributions

x <- sobolrep(model = NULL, factors = 3, N = 1000, nboot=0, nbrep=1, total=FALSE)

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
plot(x,choice=1)
plot(x,choice=2)
}


