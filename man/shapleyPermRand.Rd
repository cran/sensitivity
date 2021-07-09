\name{shapleyPermRand}
\alias{shapleyPermRand}
\alias{tell.shapleyPermRand}
\alias{print.shapleyPermRand}
\alias{plot.shapleyPermRand}
\alias{ggplot.shapleyPermRand}

\title{Estimation of Shapley effects by random permutations of inputs 
(Agorithm of Song et al, 2016), in cases of independent or dependent inputs}

\description{
  \code{shapleyPermRand} implements the Monte Carlo estimation of
  the Shapley effects (Owen, 2014) and their standard errors by randomly sampling 
  permutations of inputs (Song et al., 2016). It also estimates full first order 
  and independent total Sobol' indices (Mara et al., 2015), and their standard errors. 
  The function also allows the estimations of all these sensitivity indices in case 
  of dependent inputs.
  The total cost of   this algorithm is \eqn{Nv + m \times (d-1) \times No 
  \times Ni}{Nv + m * (d - 1) * No * Ni} model evaluations. 
}

\usage{
shapleyPermRand(model = NULL, Xall, Xset, d, Nv, m, No = 1, Ni = 3, colnames = NULL, \dots)
\method{tell}{shapleyPermRand}(x, y = NULL, return.var = NULL, \dots)
\method{print}{shapleyPermRand}(x, \dots)
\method{plot}{shapleyPermRand}(x, ylim = c(0, 1), \dots)
\method{ggplot}{shapleyPermRand}(x, ylim = c(0, 1), title = NULL, \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{Xall}{Xall(n) is a function to generate a n-sample of a d-dimensional 
  input vector (following the required joint distribution).}
  \item{Xset}{Xset(n, Sj, Sjc, xjc) is a function to generate a n-sample of 
  a d-dimensional input vector corresponding to the indices in Sj conditional 
  on the input values xjc with the index set Sjc (following the required joint 
  distribution).}
  \item{d}{number of inputs.}
  \item{Nv}{Monte Carlo sample size to estimate the output variance.}
  \item{m}{Number of randomly sampled permutations.}
  \item{No}{Outer Monte Carlo sample size to estimate the expectation 
  of the conditional variance of the model output.}
  \item{Ni}{Inner Monte Carlo sample size to estimate the 
  conditional variance of the model output.}
  \item{colnames}{Optional: A vector containing the names of the inputs.}
  \item{x}{a list of class \code{"shapleyPermRand"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{return.var}{a vector of character strings giving further
    internal variables names to store in the output object \code{x}.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{title}{a title of the plot with ggplot() function.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
}

\value{
  \code{shapleyPermRand} returns a list of class \code{"shapleyPermRand"}, containing
  all the input arguments detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments.}
  \item{y}{the response used.}
  \item{E}{the estimation of the ouput mean.}
  \item{V}{the estimation of the ouput variance.}
  \item{Shapley}{the estimations of the Shapley effects.}
  \item{SobolS}{the estimations of the full first-order Sobol' indices.}
  \item{SobolT}{the estimations of the independent total sensitivity Sobol' indices.}
  
  Users can ask more ouput variables with the argument \code{return.var}
  (for example, the list of permutations \code{perms}).
}

\details{
  This function requires R package "gtools".

  The default values No = 1 and Ni = 3 are the optimal ones obtained by the
  theoretical analysis of Song et al., 2016.
  
  The computations of the standard errors do not consider the samples to estimate 
  expectation of conditional variances. They are only made regarding the random 
  permutations and are based on the variance of the Monte carlo estimates divided 
  by m. The confidence intervals at 95\% correspond to +- 1.96 standard deviations.
}

\references{
B. Iooss and C. Prieur, 2019, \emph{Shapley effects for sensitivity analysis with 
correlated inputs: comparisons with Sobol' indices, numerical estimation and 
applications}, International Journal of Uncertainty Quantification, 9, 493--514.

T. Mara, S. Tarantola, P. Annoni, 2015, \emph{Non-parametric methods for global 
sensitivity analysis of model output with dependent inputs}, Environmental
Modeling & Software 72, 173--183.

A.B. Owen, 2014, \emph{Sobol' indices and Shapley value}, SIAM/ASA 
Journal of Uncertainty Quantification, 2, 245--251.

A.B. Owen and C. Prieur, 2016, \emph{On Shapley value for measuring importance
of dependent inputs}, SIAM/ASA Journal of Uncertainty Quantification, 5, 986--1002.

E. Song, B.L. Nelson, and J. Staum, 2016, \emph{Shapley effects for 
global sensitivity analysis: Theory and computation}, SIAM/ASA Journal 
of Uncertainty Quantification, 4, 1060--1083.

}

\author{
Bertrand Iooss, Eunhye Song, Barry L. Nelson, Jeremy Staum 
}

\seealso{
  \code{\link{shapleyPermEx}, \link{shapleyLinearGaussian}, \link{shapleySubsetMc}, \link{shapleysobol_knn}
}}

\examples{

\donttest{

##################################
# Test case : the Ishigami function
# See Iooss and Prieur (2019)

library(gtools)

d <- 3
Xall <- function(n) matrix(runif(d*n,-pi,pi),nc=d)
Xset <- function(n, Sj, Sjc, xjc) matrix(runif(n*length(Sj),-pi,pi),nc=length(Sj))

x <- shapleyPermRand(model = ishigami.fun, Xall=Xall, Xset=Xset, d=d, Nv=1e4, m=1e4, No = 1, Ni = 3)
print(x)
plot(x)

library(ggplot2)
ggplot(x)

##################################
# Test case : Linear model (3 Gaussian inputs including 2 dependent)
# See Iooss and Prieur (2019)

library(ggplot2)
library(gtools)
library(mvtnorm) # Multivariate Gaussian variables
library(condMVNorm) # Conditional multivariate Gaussian variables

modlin <- function(X) apply(X,1,sum)

d <- 3
mu <- rep(0,d)
sig <- c(1,1,2)
ro <- 0.9
Cormat <- matrix(c(1,0,0,0,1,ro,0,ro,1),d,d)
Covmat <- ( sig \%*\% t(sig) ) * Cormat

Xall <- function(n) mvtnorm::rmvnorm(n,mu,Covmat)

Xset <- function(n, Sj, Sjc, xjc){
  if (is.null(Sjc)){
    if (length(Sj) == 1){ rnorm(n,mu[Sj],sqrt(Covmat[Sj,Sj]))
    } else{ mvtnorm::rmvnorm(n,mu[Sj],Covmat[Sj,Sj])}
  } else{ condMVNorm::rcmvnorm(n, mu, Covmat, dependent.ind=Sj, given.ind=Sjc, X.given=xjc)}}

m <- 1e3 # put m)1e4 for more precised results
x <- shapleyPermRand(model = modlin, Xall=Xall, Xset=Xset, d=d, Nv=1e3, m = m, No = 1, Ni = 3)
print(x)
ggplot(x)

}

}

\keyword{design}
