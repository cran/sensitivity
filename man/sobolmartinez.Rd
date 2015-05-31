\name{sobolmartinez}
\alias{sobolmartinez}
\alias{tell.sobolmartinez}
\alias{print.sobolmartinez}
\alias{plot.sobolmartinez}

\title{Monte Carlo Estimation of Sobol' Indices (formulas of Martinez (2011))}

\description{
  \code{sobolmartinez} implements the Monte Carlo estimation of
  the Sobol' indices for both first-order and total indices using 
  correlation coefficients-based formulas, at a total cost of 
  \eqn{(p+2) \times n}{(p + 2) * n} model evaluations. 
  These are called the Martinez estimators.
}

\usage{
sobolmartinez(model = NULL, X1, X2, nboot = 0, conf = 0.95, \dots)
\method{tell}{sobolmartinez}(x, y = NULL, return.var = NULL, \dots)
\method{print}{sobolmartinez}(x, \dots)
\method{plot}{sobolmartinez}(x, ylim = c(0, 1), \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{X1}{the first random sample.}
  \item{X2}{the second random sample.}
  \item{nboot}{the number of bootstrap replicates, or zero to use theoretical 
  formulas based on confidence interfaces of correlation coefficient 
  (Martinez, 2011).}
  \item{conf}{the confidence level for bootstrap confidence intervals.}
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
  \code{sobolmartinez} returns a list of class \code{"sobolmartinez"}, 
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
  
  Users can ask more ouput variables with the argument
  \code{return.var} (for example, bootstrap outputs \code{V.boot},
  \code{S.boot} and \code{T.boot}).
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
}

\author{
Bertrand Iooss
}

\seealso{
  \code{\link{sobol}, \link{sobol2002}, \link{sobol2007}, \link{soboljansen}, 
        \link{sobolEff}, \link{sobolmara},\link{sobolMultOut}}
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

x <- sobolmartinez(model = sobol.fun, X1, X2, nboot = 0)
print(x)
plot(x)
}

\keyword{design}
