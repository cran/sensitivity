\name{delsa}
\alias{delsa}
\alias{tell.delsa}
\alias{print.delsa}
\alias{plot.delsa}

\title{Distributed Evaluation of Local Sensitivity Analysis}

\description{ \code{delsa} implements Distributed Evaluation of 
Local Sensitivity Analysis to calculate first order parameter 
sensitivity at multiple locations in parameter space. The locations 
in parameter space can either be obtained by a call to \code{\link{parameterSets}} 
or by specifying \code{X0} directly, in which case the prior variance 
of each parameter \code{varprior} also needs to be specified. Via \code{plot} 
(which uses functions of the package \code{ggplot2} and \code{reshape2}), 
the indices can be visualized.
}
\usage{
delsa(model = NULL, perturb=1.01,
  par.ranges, samples, method,
  X0, varprior,
  \dots)
  
\method{tell}{delsa}(x, y = NULL,\dots)

\method{print}{delsa}(x, \dots)

\method{plot}{delsa}(x, which=1:3, ask = dev.interactive(), \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{perturb}{Perturbation used to calculate sensitivity at each evaluation location}
  \item{par.ranges}{A named list of minimum and maximum parameter values}
  \item{samples}{Number of samples to generate. For the \code{"grid"} and
\code{"innergrid"}method, corresponds to the number of samples for
each parameter, and may be a vector.}
  \item{method}{Sampling scheme. See \code{\link{parameterSets}}}
  \item{X0}{Parameter values at which to evaluate sensitivity indices. 
  Can be used instead of specifying sampling \code{method}}
  \item{varprior}{Prior variance. If \code{X0} is specified, \code{varprior} 
  must also be specified.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
  \item{x}{a list of class \code{"delsa"} storing the state of the
    sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{which}{if a subset of the plots is required, specify a subset of the numbers 1:3}
  \item{ask}{logical; if TRUE, the user is asked before each plot, see \code{\link{par}(ask=.)}}
}

\value{
  \code{delsa} returns a list of class \code{"delsa"}, containing all
  the input arguments detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments.}
  \item{y}{a vector of model responses.}
  \item{delsafirst}{the first order indices for each location in \code{X0}}
}
\details{
\code{print} shows summary of the first order indices across parameter space.

\code{plot} shows: (1) the cumulative distribution function of first order 
sensitivity across parameter space, (2) variation of first order sensitivity 
in relation to model response, and (3) sensitivity in relation to parameter value.
}
\author{
Conversion for \code{sensitivity} package by Joseph Guillaume, 
based on original R code by Oldrich Rakovec
}
\references{
Rakovec, O., M. C. Hill, M. P. Clark, A. H. Weerts, A. J. Teuling, R. Uijlenhoet (2014), 
Distributed 
Evaluation of Local Sensitivity Analysis (DELSA), with application to hydrologic models, 
Water Resour. Res., 50, 1-18
}
\seealso{
  \code{\link{parameterSets}} which is used to generate points, \code{\link{sensitivity}} 
  for other methods in the package
}

\examples{
# Test case : the non-monotonic Sobol g-function
# (there are 8 factors, all following the uniform distribution on [0,1])

\dontrun{
library(randtoolbox)
x <- delsa(model=sobol.fun,
           par.ranges=replicate(8,c(0,1),simplify=FALSE),
           samples=100,method="sobol")

# Summary of sensitivity indices of each parameter across parameter space
print(x)

library(ggplot2)
library(reshape2)
x11()
plot(x)
}
}
\keyword{design}