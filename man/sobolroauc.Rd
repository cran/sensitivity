\name{sobolroauc}
\alias{sobolroauc}
\alias{tell.sobolroauc}
\alias{print.sobolroauc}
\alias{plot.sobolroauc}
\alias{ggplot.sobolroauc}


\title{Sobol' Indices estimation under inequality constraints}

\description{\code{sobolroauc} deals with the estimation of Sobol' sensitivity indices when there exists one or multiple sets of constrained factors. Constraints within a set are expressed as inequality constraints (simplex constraint). This function generalizes the procedure of Tissot and Prieur (2015) to estimate either all first-order indices or all closed second-order indices at a total cost of \eqn{2 \times N}{2*N} model evaluations. For closed second-order indices \eqn{N=q^{2}}{N=q^2} where \eqn{q \geq d-1}{q >= d-1} is a prime number denoting the number of levels of the orthogonal array, and where \eqn{d}{d} indicates the number of independent factors or sets of factors.}

\usage{
sobolroauc(model = NULL, factors, constraints = NULL, N, p = 1, order, 
            tail = TRUE, conf = 0.95, nboot = 0, \dots)
\method{tell}{sobolroauc}(x, y = NULL, \dots)
\method{print}{sobolroauc}(x, \dots)
\method{plot}{sobolroauc}(x, ylim = c(0,1), \dots)
\method{ggplot}{sobolroauc}(x, ylim = c(0,1), \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{factors}{an integer giving the number of factors, or a vector of character strings giving their names.}
  \item{constraints}{a list giving the sets of constrained factors (see "Details").}
  \item{N}{an integer giving the size of each replicated design (for a total of \eqn{2 \times N}{2*N} model evaluations).}
  \item{p}{an integer giving the number of model outputs.}
  \item{order}{an integer giving the order of the indices (1 or 2).}
  \item{tail}{a boolean specifying the method used to choose the number of levels of the orthogonal array (see "Warning messages").}
  \item{conf}{the confidence level for confidence intervals.}
  \item{nboot}{the number of bootstrap replicates.}
  \item{x}{a list of class \code{"sobolroauc"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
}

\value{
  \code{sobolroauc} returns a list of class \code{"sobolroauc"}, containing all
  the input arguments detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments (concatenation of two replicated designs).}
  \item{y}{the responses used.}
  \item{OA}{the orthogonal array constructed (\code{NULL} if \code{order}=1).}
  \item{V}{the estimations of Variances of the Conditional Expectations (VCE) with respect to each factor.}
  \item{S}{the estimations of the Sobol' indices.}

}

\details{

\code{constraints} list the sets of factors depending on each other through inequality constraints (see "Examples"). A same factor is not allowed to appear in multiple sets. Factors not appearing in \code{constraints} are assumed to be independent and follow each a uniform distribution on [0,1]. One Sobol' index is estimated for each independent factor or set of factors.

Missing values (i.e \code{NA} values) in the model responses are automatically handled by the function.

This function also supports multidimensional outputs (matrices in \code{y} or as output of \code{model}).
In this case, aggregated Sobol' indices are returned (\code{see sobolMultOut}).
}

\section{Warning messages}{
  \describe{
    \item{"The value entered for \code{N} is not the square of a prime number. It has been replaced by: "}{when \code{order}\eqn{=2}, the number of levels of the orthogonal array must be a prime number. If \code{N} is not a square of a prime number, then this warning message indicates that it was replaced depending on the value of \code{tail}. If \code{tail=TRUE} (resp. \code{tail=FALSE}) the new value of \code{N} is equal to the square of the prime number preceding (resp. following) the square root of \code{N}.}
    
    \item{"The value entered for \code{N} is not satisfying the constraint \eqn{N \geq (d-1)^2}{N >= (d-1)^2}. It has been replaced by: "}{when \code{order}\eqn{=2}, the following constraint must be satisfied \eqn{N \geq (d-1)^{2}}{N \ge (d-1)^2} where \eqn{d}{d} is the number of independent factors or sets of factors. This warning message indicates that \code{N} was replaced by the square of the prime number following (or equals to) \eqn{d-1}.}}

}

\references{
L. Devroye, 1986, Non-Uniform Random Variate Generation. Springer-Verlag.

J. Jacques, C. Lavergne and N. Devictor, 2006, Sensitivity Analysis in presence of model
uncertainty and correlated inputs. Reliability Engineering & System Safety, 91:1126-1134.

L. Gilquin, C. Prieur and E. Arnaud, 2015, \emph{Replication procedure for grouped Sobol' indices estimation in dependent uncertainty spaces}, Information and Inference, 4:354-379.

J.Y. Tissot and C. Prieur, 2015, \emph{A randomized orthogonal orray-based procedure for the estimation of first- and second-order Sobol' indices}, J. Statist. Comput. Simulation, 85:1358-1381.
}

\author{
Laurent Gilquin
}

\seealso{
  \code{\link{sobolroalhs}}, 
  \code{\link{sobolmara}}
}

\examples{
library(boot)
library(numbers)

# Test case: the non-monotonic Sobol g-function
# (there are 8 factors, all following the uniform distribution on [0,1])

# Suppose we have the inequality constraints: X1 <= X3 and X4 <= X6.

# first-order sensitivity indices
x <- sobolroauc(model = sobol.fun, factors = 8, constraints = list(c(1,3),c(4,6)), 
                N = 1000, order = 1, nboot=100)
print(x)
plot(x)

library(ggplot2)
ggplot(x)

# closed second-order sensitivity indices
x <- sobolroauc(model = sobol.fun, factors = 8, constraints = list(c(1,3),c(4,6)), 
                N = 1000, order = 2, nboot=100)
print(x)
ggplot(x)

}
\keyword{design}

