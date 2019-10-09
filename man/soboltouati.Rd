\name{soboltouati}
\alias{soboltouati}
\alias{tell.soboltouati}
\alias{print.soboltouati}
\alias{plot.soboltouati}
\alias{ggplot.soboltouati}

\title{Monte Carlo Estimation of Sobol' Indices (formulas of Martinez (2011) and Touati (2016))}

\description{
  \code{soboltouati} implements the Monte Carlo estimation of
  the Sobol' indices for both first-order and total indices using 
  correlation coefficients-based formulas, at a total cost of 
  \eqn{(p+2) \times n}{(p + 2) * n} model evaluations. 
  These are called the Martinez estimators. It also computes their
  confidence intervals based on asymptotic properties of empirical 
  correlation coefficients.
}

\usage{
soboltouati(model = NULL, X1, X2, conf = 0.95, \dots)
\method{tell}{soboltouati}(x, y = NULL, return.var = NULL, \dots)
\method{print}{soboltouati}(x, \dots)
\method{plot}{soboltouati}(x, ylim = c(0, 1), \dots)
\method{ggplot}{soboltouati}(x, ylim = c(0, 1), \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{X1}{the first random sample.}
  \item{X2}{the second random sample.}
  \item{conf}{the confidence level for confidence intervals, or zero to
  avoid their computation if they are not needed.}
  \item{x}{a list of class \code{"sobol"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{return.var}{a vector of character strings giving further
    internal variables names to store in  the output object \code{x}.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called}
}

\value{
  \code{soboltouati} returns a list of class \code{"soboltouati"}, 
  containing all the input arguments detailed before, plus the following 
  components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments.}
  \item{y}{the response used}
  \item{V}{the estimations of normalized variances of the Conditional 
    Expectations (VCE) with respect to each factor and also with respect 
    to the complementary set of each factor ("all but \eqn{X_i}{Xi}").}
  \item{S}{the estimations of the Sobol' first-order indices.}
  \item{T}{the estimations of the Sobol' total sensitivity indices.}
  
}

\details{
This estimator supports missing values (NA or NaN) which can occur during the
simulation of the model on the design of experiments (due to code failure) 
even if Sobol' indices are no more rigorous variance-based sensitivity 
indices if missing values are present. In this case, a warning is displayed. 
}

\references{
J-M. Martinez, 2011, \emph{Analyse de sensibilite globale par decomposition 
de la variance}, Presentation in the meeting of GdR Ondes and GdR MASCOT-NUM, 
January, 13th, 2011, Institut Henri Poincare, Paris, France.

T. Touati, 2016, Confidence intervals for Sobol' indices. 
Proceedings of the SAMO 2016 Conference, Reunion Island, France, December 2016.

T. Touati, 2017, \emph{Intervalles de confiance pour les indices de Sobol},
49emes Journees de la SFdS, Avignon, France, Juin 2017.
}

\author{
Taieb Touati, Khalid Boumhaout
}

\seealso{
  \code{\link{sobol}, \link{sobol2002}, \link{sobolSalt}, \link{sobol2007}, \link{soboljansen}, \link{sobolmartinez}}
}

\examples{
# Test case : the non-monotonic Sobol g-function

# The method of sobol requires 2 samples
# There are 8 factors, all following the uniform distribution
# on [0,1]

library(boot)
n <- 1000
X1 <- data.frame(matrix(runif(8 * n), nrow = n))
X2 <- data.frame(matrix(runif(8 * n), nrow = n))

# sensitivity analysis

x <- soboltouati(model = sobol.fun, X1, X2)
print(x)
plot(x)

library(ggplot2)
ggplot(x)
}

\keyword{design}
