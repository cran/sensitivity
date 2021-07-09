\name{lmg}
\alias{lmg}
\alias{print.lmg}
\alias{plot.lmg}

\title{LMG \eqn{R^2}{R-squared} decomposition for linear and logistic regression models}

\description{
  \code{lmg} computes the Lindeman, Merenda and Gold (LMG) indices for correlated
  input relative importance by \eqn{R^2}{R-squared} decomposition for linear and logistic
  regression models. These indices allocates  a share of \eqn{R^2}{R-squared} to each input
  based on the Shapley attribution system, in the case of dependent or correlated inputs.
}

\usage{
lmg(X, y, logistic = FALSE,  rank = FALSE, nboot = 0, 
    conf = 0.95, max.iter = 1000, parl = NULL)
\method{print}{lmg}(x, \dots)
\method{plot}{lmg}(x, ylim = c(0,1), \dots)
}

\arguments{
  \item{X}{a matrix or data frame containing the observed covariates
  (i.e., features, input variables...).}
  \item{y}{a numeric vector containing the observed outcomes (i.e.,
  dependent variable). If \code{logistic=TRUE}, can be a numeric vector
  of zeros and ones, or a logical vector, or a factor.}
  \item{logistic}{logical. If \code{TRUE}, the analysis is done via a
  logistic regression(binomial GLM).}
  \item{rank}{logical. If \code{TRUE}, the analysis is done on the
    ranks.}
  \item{nboot}{the number of bootstrap replicates for the computation
  of confidence intervals.}
  \item{conf}{the confidence level of the bootstrap confidence intervals.}
  \item{max.iter}{if \code{logistic=TRUE}, the maximum number of iterative 
  optimization steps allowed for the logistic regression. Default is \code{1000}.} 
  \item{parl}{number of cores on which to parallelize the computation. If
  \code{NULL}, then no parallelization is done.}
  \item{x}{the object returned by \code{lmg}.}
  \item{ylim}{the y-coordinate limits of the plot.}
  \item{\dots}{arguments to be passed to methods, such as graphical
    parameters (see \code{par}).}
}

\value{
  \code{lmg} returns a list of class \code{"lmg"}, containing the following
  components:

  \item{call}{the matched call.}
  \item{lmg}{the estimations of the LMG indices.}
  \item{R2s}{the estimations of the \eqn{R^2}{R-squared} for all possible sub-models.}
  \item{indices}{list of all subsets corresponding to the structure of R2s.}
  \item{w}{the Shapley weights.}
  \item{conf_int}{a matrix containing the estimations, biais and confidence
  intervals by bootstrap (if \code{nboot>0}).}
  \item{X}{the observed covariates.}
  \item{y}{the observed outcomes.}
  \item{logistic}{logical. \code{TRUE} if the analysis has been made by
  logistic regression.}
  \item{boot}{logical. \code{TRUE} if bootstrap estimates have been produced.}
  \item{nboot}{number of bootstrap replicates.}
  \item{rank}{logical. \code{TRUE} if a rank analysis has been made.}
  \item{parl}{number of chosen cores for the computation.}
  \item{conf}{level for the confidence intervals by bootstrap.}
}

\details{
  The computation is done using the subset procedure, defined in Broto, Bachoc
  and Depecker (2020), that is computing all the \eqn{R^2}{R-squared} for all possible
  sub-models first, and then affecting the Shapley weights according to the Lindeman,
  Merenda and Gold (1980) definition.

  For logistic regression (\code{logistic=TRUE}), the \eqn{R^2}{R-squared}
  value is equal to:
  \deqn{R^2 = 1-\frac{\textrm{model deviance}}{\textrm{null deviance}}}{R-squared
  = 1 - (model deviance)/(null deviance)}

  If either a logistic regression model (\code{logistic = TRUE}), or any column
  of \code{X} is categorical (i.e., of class \code{factor}), then the rank-based
  indices cannot be computed. In both those cases, \code{rank = FALSE} is forced
  by default (with a \code{warning}).

  If too many cores for the machine are passed on to the \code{parl} argument,
  the chosen number of cores is defaulted to the available cores minus one.
}

\references{

  Broto B., Bachoc F. and Depecker M. (2020) \emph{Variance Reduction for Estimation
  of Shapley Effects and Adaptation to Unknown Input Distribution.} SIAM/ASA Journal
  on Uncertainty Quantification, 8(2).
  
  D.V. Budescu (1993). \emph{Dominance analysis: A new approach to the problem of relative
  importance of predictors in multiple regression.} Psychological Bulletin, 114:542-551.

  U. Gromping (2006). \emph{Relative importance for linear regression in R: the Package
  relaimpo.}  Journal of Statistical Software, 17:1-27.

  M. Il Idrissi, V. Chabridon and B. Iooss (2021). \emph{Developments and applications
  of Shapley effects   to reliability-oriented sensitivity analysis with correlated inputs.}
  Preprint, \url{https://hal.archives-ouvertes.fr/hal-03106452}
  
  M. Il Idrissi, V. Chabridon and B. Iooss (2021). \emph{Mesures d'importance relative  
  par decompositions de la performance de modeles de regression,} Preprint, 52emes Journees 
  de Statistiques de la Societe Francaise de Statistique (SFdS), Nice, France, Juin 2021
  
  Lindeman RH, Merenda PF, Gold RZ (1980). \emph{Introduction to Bivariate and Multivariate
  Analysis.} Scott, Foresman, Glenview, IL.
}

\author{
Marouane Il Idrissi
}

\examples{
library(parallel)
library(doParallel)
library(foreach)
library(gtools)
library(boot)

library(mvtnorm)

set.seed(1234)
n <- 1000
beta<-c(1,-1,0.5)
sigma<-matrix(c(1,0,0,
                0,1,-0.8,
                0,-0.8,1),
              nrow=3,
              ncol=3)

############################
# Gaussian correlated inputs

X <-rmvnorm(n, rep(0,3), sigma)
colnames(X)<-c("X1","X2", "X3")

#############################
# Linear Model

y <- X\%*\%beta + rnorm(n,0,2)

# Without Bootstrap confidence intervals
x<-lmg(X, y)
print(x)
plot(x)

# With Boostrap confidence intervals
x<-lmg(X, y, nboot=100, conf=0.95)
print(x)
plot(x)

# Rank-based analysis
x<-lmg(X, y, rank=TRUE, nboot=100, conf=0.95)
print(x)
plot(x)

############################
# Logistic Regression
y<-as.numeric(X\%*\%beta + rnorm(n)>0)
x<-lmg(X,y, logistic = TRUE)
plot(x)
print(x)

\donttest{
# Parallel computing
x<-lmg(X,y, logistic = TRUE, parl=2)
plot(x)
print(x)
}

}

\seealso{
\code{\link{pcc}}, \code{\link{src}}, \code{\link{shapleyPermEx}}, \code{\link{shapleysobol_knn}}, \code{\link{emvd}}
}

\keyword{shapley}
