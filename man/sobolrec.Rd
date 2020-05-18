\name{sobolrec}
\alias{sobolrec}
\alias{ask.sobolrec}
\alias{tell.sobolrec}
\alias{print.sobolrec}
\alias{plot.sobolrec}


\title{Recursive estimation of Sobol' indices}

\description{\code{sobolrec} implements a recursive version of the procedure introduced by Tissot & Prieur (2015) using two replicated nested designs. This function estimates either all first-order indices or all closed second-order indices at a total cost of \eqn{2 \times N}{2*N} model evaluations where \eqn{N}{N} is the size of each replicated nested design.}

\usage{
sobolrec(model=NULL, factors, layers, order, precision, method=NULL, tail=TRUE, \dots)
\method{ask}{sobolrec}(x, index, \dots)
\method{tell}{sobolrec}(x, y = NULL, index, \dots)
\method{print}{sobolrec}(x, \dots)
\method{plot}{sobolrec}(x, ylim = c(0,1), \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{factors}{an integer giving the number of factors, or a vector of character strings giving their names.}
  \item{layers}{If \code{order=1}, a vector specifying the respective sizes of each layer (see "Details"). If \code{order=2}, an integer specifying the size of all layers.}
  \item{order}{an integer specifying which indices to estimate: \code{1} for first-order indices, \code{2} for closed second-order indices.}
  \item{precision}{a vector containing:
  \itemize{
  \item the target precision for the stopping criterion.
  \item the number of steps for the stopping criterion (must be greater than 1).}}
    \item{tail}{a boolean specifying the method used to choose the number of levels of the orthogonal array (see "Warning messages").}
  \item{method}{If \code{order=2}, a character specifying the method to construct the orthogonal arrays (see "Details"):
  \itemize{
  \item \code{"al"} for the algebraic method 
  \item \code{"ar"} for the accept-reject method
  }
  Set to \code{NULL} if \code{order=1}.}
  \item{x}{a list of class \code{"sobolrec"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{index}{an integer specifying the step of the recursion}
  \item{y}{the model response.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
}

\value{
  \code{sobolrec} returns a list of class \code{"sobolrec"}, containing all
  the input arguments detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments (row concatenation of the two replicated designs).}
  \item{y}{a list of the response used at each step.}
  \item{V}{a list of the model variance estimated at each step.}
  \item{S}{a list of the Sobol' indices estimated at each step.}
  \item{steps}{the number of steps performed.}
  \item{N}{the size of each replicated nested design.}
}

\details{
For first-order indices, \code{layers} is a vector:
\deqn{\left(s_1, ...,s_m \right)}{(s_1, ..., s_m)}
specifying the number \eqn{m}{m} of layers of the nested design whose respective size are given by: 
\deqn{\prod_{i=1}^{k-1} s_i, \ k=2, ...,m+1}{s_1 * s_2 * ... * s_{k-1}, for k=2, ...,m+1.}
For closed second-order indices, \code{layers} directly specifies the size of all layers.

For each Sobol' index \eqn{S}{S} the stopping criterion writes:
\deqn{\mid S_{l-1}-S_{l} \mid < \epsilon}{| S_{l-1}-S_{l} | < \epsilon}
This criterion is tested for the last \eqn{l_0}{l_0} steps (including the current one). \eqn{\epsilon}{\epsilon} and \eqn{l_0}{l_0} are respectively the target precision and the number of steps of the stopping criterion specified in \code{precision}.

\code{sobolrec} uses either an algebraic or an accept-rejet \code{method} to construct the orthogonal arrays for the estimation of closed second-order indices. The algebraic method is less precise than the accept-reject method but offers more steps when the number of \code{factors} is small.


\code{sobolrec} automatically assigns a uniform distribution on [0,1] to each input. Transformations of distributions (between U[0,1] and the wanted distribution) have to be performed before the call to tell().
}

\section{Warning messages}{
  \describe{
    \item{"The value entered for \code{layers} is not the square of a prime number. It has been replaced by: "}{When \code{order=2}, the value of \code{layers} must be the square of a prime power number.  This warning message indicates that it was not the case and the value has been replaced depending on \code{tail}. If \code{tail=TRUE} (resp. \code{tail=FALSE}) the new value of \code{layers} is equal to the square of the prime number preceding (resp. following) the square root of \code{layers}.}
    
    \item{"The value entered for \code{layers} is not satisfying the constraint. It has been replaced by: "}{the value \eqn{N} for \code{layers} must satisfied the constraint \eqn{N \geq (d-1)^{2}}{N \ge (d-1)^2} where \eqn{d}{d} is the number of factors. This warning message indicates that \code{N} was replaced by the square of the prime number following (or equals to) \eqn{d-1}.}}

}

\references{
A.S. Hedayat, N.J.A. Sloane and J. Stufken, 1999, \emph{Orthogonal Arrays: Theory and Applications}, Springer Series in Statistics.

L. Gilquin, E. Arnaud, H. Monod and C. Prieur, 2016, \emph{Recursive estimation procedure of Sobol' indices based on replicated designs}, preprint available at: \url{https://hal.inria.fr/hal-01291769}.
}


\examples{
	\dontrun{
# Test case: the non-monotonic Sobol g-function

# The method of sobol requires 2 samples
# (there are 8 factors, all following the uniform distribution on [0,1])

# first-order indices estimation
x <- sobolrec(model = sobol.fun, factors = 8, layers=rep(2,each=15), order=1,
              precision = c(5*10^(-2),2), method=NULL, tail=TRUE)
print(x)

# closed second-order indices estimation
x <- sobolrec(model = sobol.fun, factors = 8, layers=11^2, order=2,
              precision = c(10^(-2),3), method="al", tail=TRUE)
print(x)


# Test case: dealing with external model

x <- sobolrec(model = NULL, factors = 8, layers=rep(2,each=15), order=1,
              precision = c(5*10^(-2),2), method=NULL, tail=TRUE)
toy <- sobol.fun
k <- 1
stop_crit <- FALSE
while(!(stop_crit) & (k<length(x$layers))){
  ask(x, index=k)
  y <- toy(x$block)
  tell(x, y, index=k)
  stop_crit <- x$stop_crit
  k <- k+1
}
print(x)
}
}

