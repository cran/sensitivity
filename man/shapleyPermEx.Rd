\name{shapleyPermEx}
\alias{shapleyPermEx}
\alias{tell.shapleyPermEx}
\alias{print.shapleyPermEx}
\alias{plot.shapleyPermEx}
\alias{ggplot.shapleyPermEx}

\title{Estimation of Shapley effects by examining all permutations of inputs 
(Agorithm of Song et al, 2016), in cases of independent or dependent inputs}

\description{
  \code{shapleyPermEx} implements the Monte Carlo estimation of
  the Shapley effects (Owen, 2014) and their standard errors by examining all
  permutations of inputs (Song et al., 2016; Iooss and Prieur, 2019). It also 
  estimates full first order and independent total Sobol' indices 
  (Mara et al., 2015). The function also allows the estimations of all these 
  sensitivity indices in case of dependent inputs. The total cost of this 
  algorithm is \eqn{Nv + d! \times (d-1) \times No \times Ni}{Nv + d! * (d - 1) * No * Ni} 
  model evaluations. 
}

\usage{
shapleyPermEx(model = NULL, Xall, Xset, d, Nv, No, Ni = 3, colnames = NULL, \dots)
\method{tell}{shapleyPermEx}(x, y = NULL, return.var = NULL, \dots)
\method{print}{shapleyPermEx}(x, \dots)
\method{plot}{shapleyPermEx}(x, ylim = c(0, 1), \dots)
\method{ggplot}{shapleyPermEx}(data, mapping = aes(), ylim = c(0, 1), title = NULL,
                 \dots, environment = parent.frame())
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
  \item{No}{Outer Monte Carlo sample size to estimate the expectation 
  of the conditional variance of the model output.}
  \item{Ni}{Inner Monte Carlo sample size to estimate the 
  conditional variance of the model output.}
  \item{colnames}{Optional: A vector containing the names of the inputs.}
  \item{x}{a list of class \code{"shapleyPermEx"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{data}{a list of class \code{"shapleyPermEx"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{return.var}{a vector of character strings giving further
    internal variables names to store in the output object \code{x}.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{title}{a title of the plot with ggplot() function.}
  \item{mapping}{Default list of aesthetic mappings to use for plot. If not specified, 
    must be supplied in each layer added to the plot.}
  \item{environment}{[Deprecated] Used prior to tidy evaluation.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
}

\value{
  \code{shapleyPermEx} returns a list of class \code{"shapleyPermEx"}, containing
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

  The default values Ni = 3 is the optimal one obtained by the theoretical
  analysis of Song et al., 2016.
  
  The computations of the standard errors (and then the confidence intervals)
  come from Iooss and prieur (2019). Based on the outer Monte carlo loop 
  (calculation of expectation of conditional variance), the variance of the 
  Monte carlo estimate is divided by No. The standard error is then averaged 
  over the exact permutation loop. The confidence intervals at 95\% correspond 
  to +- 1.96 standard deviations.
}

\references{
B. Iooss and C. Prieur, 2019, \emph{Shapley effects for sensitivity analysis with 
correlated inputs: comparisons with Sobol' indices, numerical estimation and 
applications}, International Journal for Uncertainty Quantification, 9, 493--514.

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
  \code{\link{shapleyPermRand}, \link{shapleyLinearGaussian},  \link{shapleySubsetMc}, \link{shapleysobol_knn}}, \code{\link{lmg}}
}

\examples{

\donttest{

##################################
# Test case : the Ishigami function (3 uniform independent inputs)
# See Iooss and Prieur (2019)

library(gtools)

d <- 3
Xall <- function(n) matrix(runif(d*n,-pi,pi),nc=d)
Xset <- function(n, Sj, Sjc, xjc) matrix(runif(n*length(Sj),-pi,pi),nc=length(Sj))

x <- shapleyPermEx(model = ishigami.fun, Xall=Xall, Xset=Xset, d=d, Nv=1e4, No = 1e3, Ni = 3)
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
  } else{ condMVNorm::rcmvnorm(n, mu, Covmat, dependent.ind=Sj, given.ind=Sjc, 
                                X.given=xjc)}}

x <- shapleyPermEx(model = modlin, Xall=Xall, Xset=Xset, d=d, Nv=1e4, 
                    No = 1e3, Ni = 3)
print(x)
ggplot(x)
}

}
\keyword{design}
