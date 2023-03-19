\name{discrepancyCriteria_cplus}
\alias{discrepancyCriteria_cplus}

\title{Discrepancy measure}
\description{Compute discrepancy criteria. This function uses a C++ implementation of the function \code{discrepancyCriteria} from package \pkg{DiceDesign}.}
\usage{discrepancyCriteria_cplus(design,type='all')}

\arguments{
  \item{design}{a matrix corresponding to the design of experiments.
  The discrepancy criteria are computed for a design in the unit cube [0,1]\eqn{^d}.
  If this condition is not satisfied the design is automatically rescaled.}
  \item{type}{type of discrepancies (single value or vector) to be computed:
	\tabular{ll}{
		\code{'all'} \tab all type of discrepancies (default) \cr
  	\code{'C2'} \tab centered L2-discrepancy \cr
  	\code{'L2'} \tab L2-discrepancy \cr
  	\code{'L2star'} \tab L2star-discrepancy \cr
		\code{'M2'} \tab modified L2-discrepancy \cr
		\code{'S2'} \tab symmetric L2-discrepancy \cr
		\code{'W2'} \tab wrap-around L2-discrepancy \cr
	}
  }
}

\details{The discrepancy measures how far a given distribution of points deviates
from a perfectly uniform one. Different discrepancies are available.
For example, if we denote by \eqn{Vol(J)} the volume of a subset \eqn{J} of \eqn{[0; 1]^d} and \eqn{A(X; J)} the number of points of \eqn{X} falling in \eqn{J}, the \eqn{L2} discrepancy is:
\deqn{D_{L2} (X) = \left[ \int_{[0,1]^{2d}}{} \left( \frac{A(X,J_{a,b})}{n} - Vol (J_{a,b}) \right)^{2} da db \right]^{1/2}}{DL2 (X)^2 = \int_{[0,1]^d} [(A(X,J_{a,b})/n - Vol(J_{a,b})]^2 da db}
where \eqn{a = (a_{1}; ... ; a_{d})'}{a = (a1; ... ; ad)'}, \eqn{b = (b_{1};...; b_{d})'}{b = (b1;...; bd)'} and \eqn{J_{a,b} =
[a_{1}; b_{1}) \times ... \times [a_{d};b_{d})}{J_{a,b} =
[a1; b1) X ... X [ad;bd)}. The other L2-discrepancies are defined according to the same principle with different form from the subset \eqn{J}.
Among all the possibilities, discrepancyCriteria_cplus implements only the L2 discrepancies because it can be expressed analytically even for high dimension.


Centered L2-discrepancy is computed using the analytical expression done by Hickernell (1998). The user will refer to Pleming and Manteufel (2005) to have more details about the wrap around discrepancy.
}

\value{A list containing the L2-discrepancies of the \code{design}.}

\references{
Fang K.T, Li R. and Sudjianto A. (2006) Design and Modeling for
	Computer Experiments, \emph{Chapman & Hall}.

Hickernell F.J. (1998) A generalized discrepancy and quadrature error bound.
	\emph{Mathematics of Computation}, \bold{67}, 299-322.

Pleming J.B. and Manteufel R.D. (2005) \emph{Replicated Latin Hypercube Sampling},
46th Structures, Structural Dynamics & Materials Conference, 16-21 April 2005, Austin
(Texas) -- AIAA 2005-1819.
}

\author{Laurent Gilquin}

\seealso{The distance criterion provided by \code{\link{maximin_cplus}}}

\examples{
dimension <- 2
n <- 40
X <- matrix(runif(n*dimension),n,dimension)
discrepancyCriteria_cplus(X)
}

\keyword{ design }

