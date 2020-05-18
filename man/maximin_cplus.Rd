\name{maximin_cplus}
\alias{maximin_cplus}
\title{Maximin criterion}
\description{Compute the \code{maximin} criterion (also called mindist). This function uses a C++ implementation of the function mindist from package \pkg{DiceDesign}.}
\usage{maximin_cplus(design)}

\arguments{
  \item{design}{a matrix representing the design of experiments in the unit cube [0,1]\eqn{^d}. If this last condition is not fulfilled, a transformation into [0,1]\eqn{^{d}} is applied before the computation of the criteria.}
}

\details{The maximin criterion is defined by:
\deqn{maximin= \min_{x_{i}\in X} \left( \gamma_{i} \right)}{maximin = min (g_1, ... g_n)}
where \eqn{\gamma_{i}}{g_i} is the minimal distance between the point \eqn{x_{i}}{x_i}
and the other points \eqn{x_{k}}{x_k} of the \code{design}.

A higher value corresponds to a more regular scaterring of design points.
}

\value{A real number equal to the value of the maximin criterion for the \code{design}.}

\references{Gunzburer M., Burkdart J. (2004) \emph{Uniformity measures for point
samples in hypercubes}
\url{https://people.sc.fsu.edu/~jburkardt/}.

Jonshon M.E., Moore L.M. and Ylvisaker D. (1990) \emph{Minmax and maximin distance
designs}, J. of Statis. Planning and Inference, 26, 131-148.

Chen V.C.P., Tsui K.L., Barton R.R. and Allen J.K. (2003) \emph{A review of design
and modeling in computer experiments}, Handbook of Statistics, 22, 231-261.
}

\author{Laurent Gilquin}

\seealso{discrepancy measures provided by \code{\link{discrepancyCriteria_cplus}}.}

\examples{
dimension <- 2
n <- 40
X <- matrix(runif(n*dimension),n,dimension)
maximin_cplus(X)
}
