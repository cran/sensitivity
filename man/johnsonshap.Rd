\name{johnsonshap}
\alias{johnsonshap}
\alias{print.johnsonshap}
\alias{plot.johnsonshap}
\alias{ggplot.johnsonshap}

\title{Johnson-Shapley indices}

\description{
  \code{johnsonshap} computes the Johnson-Shapley indices for correlated input 
  relative importance. These indices allocate a share of the output variance to 
  each input based on the relative weight allocation system, 
  in the case of dependent or correlated inputs.
}

\usage{
johnsonshap(model = NULL, X1, N, nboot = 0, conf = 0.95)
\method{print}{johnsonshap}(x, \dots)
\method{plot}{johnsonshap}(x, ylim = c(0,1), \dots)
\method{ggplot}{johnsonshap}(data,  mapping = aes(), ylim = c(0, 1), \dots, 
                environment = parent.frame())
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{X1}{a data frame (or object coercible by \code{as.data.frame})
    containing a design of experiments (model input variables).}
  \item{N}{an integer giving the size of each replicated design for
  the Sobol' indices computations via the sobolrep() fct.}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level of the bootstrap confidence intervals.}
  \item{x}{the object returned by \code{johnsonshap}.}
  \item{data}{the object returned by \code{johnsonshap}.}
  \item{ylim}{the y-coordinate limits of the plot.}
  \item{mapping}{Default list of aesthetic mappings to use for plot. If not 
    specified, must be supplied in each layer added to the plot.}
  \item{environment}{[Deprecated] Used prior to tidy evaluation.}
  \item{\dots}{arguments to be passed to methods, such as graphical
    parameters (see \code{par}).}
}

\value{
  \code{johnsonshap} returns a list of class \code{"johnsonshap"}, containing 
  all the input arguments detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a matrix containing the design of experiments.}
  \item{sobrepZ}{the Sobol' indices of the transformed inputs (independent)}
  \item{Wstar}{the standardized weight matrix.}
  \item{johnsonshap}{a data frame containing the estimations of the 
    Johnson-Shapley indices, bias and confidence intervals.}
}

\details{
  X1 is not used to run the model but just to perform the SVD; the model is run 
  on a specific   design which is internally generated.
  
  By using bootstrap, values in the columns 'bias' and 'std. error' are 
  arbitrarily put at 0 because of impossible computations; values in columns 
  'original', 'min c.i.' and 'max c.i.' are correctly computed.
}

\references{

B. Iooss and L. Clouvel, \emph{Une methode d'approximation des effets 
de Shapley en grande dimension}, 54emes  Journees de Statistique, 
Bruxelles, Belgique, July 3-7, 2023

}

\author{
Bertrand Iooss
}

\examples{

library(ggplot2)
library(boot)

#####################################################
# Test case: the non-monotonic Sobol g-function (with independent inputs)
n <- 1000
X <- data.frame(matrix(runif(8 * n), nrow = n))
x <- johnsonshap(model = sobol.fun, X1 = X, N = n)
print(x)
plot(x)
ggplot(x)

\donttest{
#############################################
# 3D analytical toy functions described in Iooss & Clouvel (2023)

library(mvtnorm)

Xall <- function(n) mvtnorm::rmvnorm(n,mu,Covmat)
# 2 correlated inputs
Cov3d2 <- function(rho){ # correl (X1,X2)
  Cormat <- matrix(c(1,rho,0,rho,1,0,0,0,1),3,3)
  return( ( sig \%*\% t(sig) ) * Cormat)
}
mu3d <- c(1,0,0) ; sig3d <- c(0.25,1,1)
d <- 3 ; mu <- mu3d ; sig <- sig3d ; Covm <- Cov3d2
Xvec <- c("X1","X2","X3")

n <- 1e4    # initial sample size
N <- 1e4    # cost to estimate indices 
rho <- 0.9  # correlation coef for dependent inputs' case

################
# Linear model + a strong 2nd order interaction

toy3d <- function(x) return(x[,1]*(1+x[,1]*(cos(x[,2]+x[,3])^2))) 
# interaction X2X3
toy <- toy3d 

# Independent case

Covmat <- Covm(0)
X <- as.data.frame(Xall(n))
Y <- toy(X)
joh <- johnson(X, Y, nboot=100)
print(joh)
johshap <- johnsonshap(model = toy, X1 = X, N = N, nboot=100)
print(johshap)
ggplot(johshap)

# Dependent case

Covmat <- Covm(rho)
Xdep <- as.data.frame(Xall(n))
Ydep <- toy(Xdep)
joh <- johnson(Xdep, Ydep, nboot=0)
print(joh)
johshap <- johnsonshap(model = toy, X1 = Xdep, N = N, nboot=100)
print(johshap)
ggplot(johshap)

################
# Strongly non-inear model + a strong 2nd order interaction

toy3dNL <- function(x) return(sin(x[,1]*pi/2)*(1+x[,1]*(cos(x[,2]+x[,3])^2))) 
# non linearity in X1
toy <- toy3dNL

# Independent case

Covmat <- Covm(0)
X <- as.data.frame(Xall(n))
Y <- toy(X)
joh <- johnson(X, Y, nboot=100)
print(joh)
johshap <- johnsonshap(model = toy, X1 = X, N = N, nboot=100)
print(johshap)
ggplot(johshap)

# Dependent case

Covmat <- Covm(rho)
Xdep <- as.data.frame(Xall(n))
Ydep <- toy(Xdep)
joh <- johnson(Xdep, Ydep, nboot=0)
print(joh)
johshap <- johnsonshap(model = NULL, X1 = Xdep, N = N, nboot=100)
y <- toy(johshap$X)
tell(johshap, y)
print(johshap)
ggplot(johshap)

}
}

\seealso{
\code{\link{johnson}, \link{shapleysobol_knn}}
}

\keyword{regression}
