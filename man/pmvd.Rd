\name{pmvd}
\alias{pmvd}
\alias{print.pmvd}
\alias{plot.pmvd}

\title{Proportional Marginal Variance Decomposition indices
      for linear and logistic models}

\description{
  \code{pmvd} computes the PMVD indices derived from Feldman (2005) applied to
  the explained variance (\eqn{R^2}{R-squared}) as a performance metric. 
  They allow for relative importance indices by \eqn{R^2}{R-squared} decomposition 
  for linear and logistic regression models. These indices allocate a share of
  \eqn{R^2}{R-squared} to each input based on a Proportional attribution system,
  allowing for covariates with null regression coefficients to have indices
  equal to 0, despite their potential dependence with other covariates (Exclusion
  principle).
}

\usage{
pmvd(X, y, logistic = FALSE, tol = NULL, rank = FALSE, nboot = 0, 
    conf = 0.95, max.iter = 1000, parl = NULL)
\method{print}{pmvd}(x, \dots)
\method{plot}{pmvd}(x, ylim = c(0,1), \dots)
}

\arguments{
  \item{X}{a matrix or data frame containing the observed covariates
  (i.e., features, input variables...).}
  \item{y}{a numeric vector containing the observed outcomes (i.e.,
  dependent variable). If \code{logistic=TRUE}, can be a numeric vector
  of zeros and ones, or a logical vector, or a factor.}
  \item{logistic}{logical. If \code{TRUE}, the analysis is done via a
  logistic regression(binomial GLM).}
  \item{tol}{covariates with absolute marginal contributions less or equal to 
    \code{tol} are omitted. By default, if \code{tol=NULL}, only covariates with no 
    marginal contribution are omitted.}
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
  \code{pmvd} returns a list of class \code{"pmvd"}, containing the following
  components:

  \item{call}{the matched call.}
  \item{pmvd}{a data frame containing the estimations of the PMVD indices.}
  \item{R2s}{the estimations of the \eqn{R^2}{R-squared} for all possible sub-models.}
  \item{indices}{list of all subsets corresponding to the structure of R2s.}
  \item{P}{the values of \eqn{P(.)} of all subsets for recursive computing. Equal to
  \code{NULL} if bootstrap estimates are made.}
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

  The computation of the PMVD is done using the recursive method defined in
  Feldman (2005), but using the subset procedure defined in Broto, Bachoc
  and Depecker (2020), that is computing all the \eqn{R^2}{R-squared} for all
  possible sub-models first, and then computing \eqn{P(.)} recursively for all
  subsets of covariates. See Il Idrissi et al. (2021).

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
  
  Spurious covariates are defined by the \code{tol} argument. If \code{null},
  then covariates with:
  \deqn{w(\{i\}) = 0}{w(\{i\})=0}
  are omitted, and their \code{pmvd} index is set to zero. In other cases, the 
  spurious covariates are detected by:
  \deqn{|w(\{i\})| \leq \textrm{tol}}{|w(\{i\})| <= tol}

}

\references{

  Broto B., Bachoc F. and Depecker M. (2020) \emph{Variance Reduction for Estimation
  of Shapley Effects and Adaptation to Unknown Input Distribution.} SIAM/ASA Journal
  on Uncertainty Quantification, 8(2).
  
  D.V. Budescu (1993). \emph{Dominance analysis: A new approach to the problem of relative
  importance of predictors in multiple regression.} Psychological Bulletin, 114:542-551.

  L. Clouvel, B. Iooss, V. Chabridon, M. Il Idrissi and F. Robin, 2024,
  \emph{An overview of variance-based importance measures in the linear regression context: 
  comparative analyses and numerical tests}, Preprint.
  \url{https://hal.science/hal-04102053}
  
  Feldman, B. (2005) \emph{Relative Importance and Value} SSRN Electronic Journal.

  U. Gromping (2006). \emph{Relative importance for linear regression in R: the Package
  relaimpo.}  Journal of Statistical Software, 17:1-27.
  
  M. Il Idrissi, V. Chabridon and B. Iooss (2021). \emph{Mesures d'importance relative  
  par decompositions de la performance de modeles de regression,} Actes des 52emes Journees 
  de Statistiques de la Societe Francaise de Statistique (SFdS), pp 497-502,
  Nice, France, Juin 2021

  B. Iooss, V. Chabridon and V. Thouvenot, \emph{Variance-based importance 
  measures for machine learning model interpretability}, Congres lambda-mu23,
  Saclay, France, 10-13 octobre 2022
  \url{https://hal.science/hal-03741384}
  
}

\author{
Marouane Il Idrissi
}

\examples{
library(parallel)
library(gtools)
library(boot)

library(mvtnorm)

set.seed(1234)
n <- 100
beta<-c(1,-2,3)
sigma<-matrix(c(1,0,0,
                0,1,-0.8,
                0,-0.8,1),
              nrow=3,
              ncol=3)

############################
# Gaussian correlated inputs

X <-rmvnorm(n, rep(0,3), sigma)

#############################
# Linear Model

y <- X\%*\%beta + rnorm(n)

# Without Bootstrap confidence intervals
x<-pmvd(X, y)
print(x)
plot(x)

# With Boostrap confidence intervals
x<-pmvd(X, y, nboot=100, conf=0.95)
print(x)
plot(x)

# Rank-based analysis
x<-pmvd(X, y, rank=TRUE, nboot=100, conf=0.95)
print(x)
plot(x)

############################
# Logistic Regression
y<-as.numeric(X\%*\%beta + rnorm(n)>0)
x<-pmvd(X,y, logistic = TRUE)
plot(x)
print(x)

# Parallel computing
#x<-pmvd(X,y, logistic = TRUE, parl=2)
#plot(x)
#print(x)

}

\seealso{
\code{\link{pcc}}, \code{\link{src}}, \code{\link{lmg}}, \code{\link{pme_knn}}
}

\keyword{shapley}
