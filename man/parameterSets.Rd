\name{parameterSets}
\alias{parameterSets}
\title{
Generate parameter sets
}
\description{
Generate parameter sets from given ranges, with chosen sampling scheme
}
\usage{
parameterSets(par.ranges, samples, method = c("sobol", "innergrid", "grid"))
}
\arguments{
  \item{par.ranges}{
A named list of minimum and maximum parameter values
}
  \item{samples}{
Number of samples to generate. For the \code{"grid"} and \code{"innergrid"}
method, may be a vector of number of samples for each parameter.
}
  \item{method}{
the sampling scheme; see Details
}
}
\details{
Method \code{"sobol"} generates uniformly distributed Sobol low discrepancy numbers, 
using the sobol function in the randtoolbox package.

Method \code{"grid"} generates a grid within the parameter ranges, including its extremes, 
with number of points determined by \code{samples}

Method \code{"innergrid"} generates a grid within the parameter ranges, with edges 
of the grid offset from the extremes. The offset is calculated as half 
of the resolution of the grid \code{diff(par.ranges)/samples/2}.

}
\value{
the result is a \code{matrix}, with named columns for each parameter in \code{par.ranges}. 
Each row represents one parameter set.
}
\author{
Joseph Guillaume, based on similar function by Felix Andrews
}
\seealso{
\code{\link{delsa}}, which uses this function
}
\examples{

X.grid <- parameterSets(par.ranges=list(V1=c(1,1000),V2=c(1,4)),
                          samples=c(10,10),method="grid")
plot(X.grid)

X.innergrid<-parameterSets(par.ranges=list(V1=c(1,1000),V2=c(1,4)),
                          samples=c(10,10),method="innergrid")
points(X.innergrid,col="red")

\dontrun{
library(randtoolbox)
X.sobol<-parameterSets(par.ranges=list(V1=c(1,1000),V2=c(1,4)),
                           samples=100,method="sobol")
plot(X.sobol)
}
}
\keyword{ utilities }
