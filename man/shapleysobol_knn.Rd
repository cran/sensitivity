\name{shapleysobol_knn}
\alias{shapleysobol_knn}
\alias{tell.shapleysobol_knn}
\alias{extract.shapleysobol_knn}
\alias{print.shapleysobol_knn}
\alias{plot.shapleysobol_knn}
\alias{ggplot.shapleysobol_knn}
\alias{print.sobol_knn}
\alias{plot.sobol_knn}



\title{Data given Shapley effects estimation via nearest-neighbors procedure}

\description{
  \code{shapleysobol_knn} implements the estimation of several sensitivity indices using 
  only N model evaluations via ranking (following Gamboa et al. (2020) and Chatterjee (2019)) 
  or nearest neighbour search (Broto et al. (2020) and Azadkia & Chatterjee (2020)). 
  Parallelized computations are possible to accelerate the estimation process.
  It can be used with categorical inputs (which are transformed with one-hot encoding), 
  dependent inputs and multiple outputs. Sensitivity indices of any group of inputs can be computed,
  which means that in particular (full) first-order, (independent) total Sobol indices 
  and Shapley effects are accessible. For large sample sizes, the nearest neightbour algorithm 
  can be significantly accelerated   by using approximate nearest neighbour search. 
  It is also possible to estimate Shapley effects   with the random permutation approach of 
  Castro et al.(2009), where all the terms are obtained with ranking or nearest neighbours.
}

\usage{
shapleysobol_knn(model=NULL, X, method = "knn", n.knn = 2, n.limit = 2000, U = NULL,
          n.perm = NULL, noise = F, rescale = F, nboot = NULL, boot.level = 0.8, 
          conf=0.95, parl=NULL, \dots)
\method{tell}{shapleysobol_knn}(x, y, \dots)
\method{extract}{shapleysobol_knn}(x, \dots)
\method{print}{shapleysobol_knn}(x, \dots)
\method{plot}{shapleysobol_knn}(x, ylim = c(0,1), \dots)
\method{ggplot}{shapleysobol_knn}(data, mapping = aes(), ylim = c(0, 1), \dots, environment
                 = parent.frame())
\method{print}{sobol_knn}(x, \dots)
\method{plot}{sobol_knn}(x, ylim = c(0,1), \dots)
}

\arguments{
  \item{model}{a function defining the model to analyze, taking X as an argument.}
  \item{X}{a matrix or data frame containing the observed inputs.}
  \item{method}{the algorithm to be used for estimation, either "rank" or "knn",
  see details. Default is \code{method="knn"}.}
  \item{n.knn}{the number of nearest neighbours used for estimation.}
  \item{n.limit}{sample size limit above which approximate nearest neighbour
  search is activated.}
  \item{U}{an integer equal to 0 (total Sobol indices) or 1 (first-order Sobol indices) 
    or a list of vector indices defining the subsets of inputs whose sensitivity indices 
    must be computed or a matrix of 0s and 1s where each row encodes a subset of inputs 
    whose sensitivity indices must be computed (see examples). Default value is \code{NULL},
    meaning that Shapley values are returned (see details).}
  \item{n.perm}{an integer, indicating the number of random permutations used 
  for the Shapley effects' estimation. Default is \code{n.perm=NULL}, indicating
  that all possible permutations are used.}
  \item{noise}{a logical which is TRUE if the model or the output sample is 
  noisy. See details.}
  \item{rescale}{a logical indicating if continuous inputs must be rescaled before distance computations.
    If TRUE, continuous inputs are first whitened with the ZCA-cor whitening procedure 
    (cf. whiten() function in package \code{whitening}). If the inputs are independent, 
    this first step will have a very limited impact. Then, the resulting whitened inputs 
    are individually modified via a copula transform such that each input has the same scale.}
  \item{nboot}{the number of bootstrap resamples for the bootstrap estimate of 
  confidence intervals. See details.}
  \item{boot.level}{a numeric between 0 and 1 for the proportion of the 
  bootstrap sample size.}
  \item{conf}{the confidence level of the bootstrap confidence intervals.}
  \item{parl}{number of cores on which to parallelize the computation. If
  \code{NULL}, then no parallelization is done.}
  \item{x}{the object returned by \code{shapleysobol_knn}.}
  \item{data}{the object returned by \code{shapleysobol_knn}.}
  \item{y}{a numeric univariate vector containing the observed outputs.}
  \item{ylim}{the y-coordinate limits for plotting.}
  \item{mapping}{Default list of aesthetic mappings to use for plot. If not specified, 
    must be supplied in each layer added to the plot.}
  \item{environment}{[Deprecated] Used prior to tidy evaluation.}
  \item{\dots}{additional arguments to be passed to \code{model}, or to the 
  methods, such as graphical parameters (see \code{par}).}
}

\value{
  \code{shapleysobol_knn} returns a list of class \code{"shapleysobol_knn"} if \code{U=NULL}, 
  containing the following components:

  \item{call}{the matched call.}
  \item{Shap}{the estimations of the Shapley effect indices.}
  \item{VE}{the estimations of the closed Sobol' indices for all possible sub-models.}
  \item{indices}{list of all subsets corresponding to the structure of VE.}
  \item{method}{which estimation method has been used.}
  \item{n.perm}{number of random permutations.}
  \item{w}{the Shapley weights.}
  \item{conf_int}{a matrix containing the estimations, biais and confidence
  intervals by bootstrap (if \code{nboot>0}).}
  \item{X}{the observed covariates.}
  \item{y}{the observed outcomes.}
  \item{n.knn}{value of the \code{n.knn} argument.}
  \item{n.limit}{value of the \code{n.limit} argument.}
  \item{U}{value of the \code{U} argument.}
  \item{rescale}{wheter the design matrix has been rescaled.}
  \item{n.limit}{maximum number of sample before nearest-neighbor approximation.}
  \item{boot.level}{value of the \code{boot.level} argument.}
  \item{noise}{wheter the Shapley values must sum up to one or not.}
  \item{boot}{logical, wheter bootstrap confidence interval estimates have 
  been performed.}
  \item{nboot}{value of the \code{nboot} argument.}
  \item{parl}{value of the \code{parl} argument.}
  \item{conf}{value of the \code{conf} argument.}
  
  \code{shapleysobol_knn} returns a list of class \code{"sobol_knn"} if \code{U}, 
  is specified, containing the following components:

  \item{call}{the matched call.}
  \item{Sobol}{the estimations of the Sobol' indices.}
  \item{indices}{list of all subsets corresponding to the structure of VE.}
  \item{method}{which estimation method has been used.}
  \item{conf_int}{a matrix containing the estimations, biais and confidence
  intervals by bootstrap (if \code{nboot>0}).}
  \item{X}{the observed covariates.}
  \item{y}{the observed outcomes.}
  \item{U}{value of the \code{U} argument.}
  \item{n.knn}{value of the \code{n.knn} argument.}
  \item{rescale}{wheter the design matrix has been rescaled.}
  \item{n.limit}{value of the \code{n.limit} argument.}
  \item{boot.level}{value of the \code{boot.level} argument.}
  \item{boot}{logical, wheter bootstrap confidence interval estimates have 
  been performed.}
  \item{nboot}{value of the \code{nboot} argument.}
  \item{parl}{value of the \code{parl} argument.}
  \item{conf}{value of the \code{conf} argument.}
}

\details{
  For \code{method="rank"}, the estimator is defined in Gamboa et al. (2020) 
  following Chatterjee (2019). For first-order indices it is based on an input 
  ranking (same algorithm as in \code{sobolrank}) while for higher orders, 
  it uses an approximate heuristic solution of the traveling salesman problem 
  applied to the input sample distances (cf. TSP() function in package 
  \code{TSP}).  For \code{method="knn"}, ranking and TSP are replaced by a 
  nearest neighbour search   as proposed in Broto et al. (2020) and in Azadkia 
  & Chatterjee (2020) for a similar coefficient. 

  The computation is done using the subset procedure, defined in Broto, Bachoc
  and Depecker (2020), that is computing all the Sobol' closed indices for all 
  possible sub-models first, and then affecting the Shapley weights.

  It is the same algorithm as \code{sobolshap_knn} with \code{method = "knn"}
  with a slight computational improvement (the search for weight affectations is
  done on much smaller matrices, stored in a list indexed by their order), and 
  ability to perform parallel computation   and boostrap confidence interval 
  estimates.
  
  Since boostrap creates ties which are not accounted for in the algorithm, 
  confidence intervals are obtained by sampling without replacement with a 
  proportion of the total sample size \code{boot.level}, drawn uniformly.
  
  If the outputs are noisy, the argument \code{noise} can be used: it only has 
  an impact on the estimation of one specific sensitivity index, namely 
  \eqn{Var(E(Y|X1,\ldots,Xp))/Var(Y)}. If there is no noise this index is equal 
  to 1, while in the presence of noise it must be estimated.

  The distance used for subsets with mixed inputs (continuous and categorical) 
  is the Euclidean distance, thanks to a one-hot encoding of categorical inputs.
  
  If too many cores for the machine are passed on to the \code{parl} argument,
  the chosen number of cores is defaulted to the available cores minus one.
  
  If argument \code{U} is specified, only the estimated first-order or total 
  Sobol' indices are returned, or the estimated closed Sobol' indices for the 
  selected subsets. The Shapley effects are not computed, and thus, not returned.
  
  The \code{extract} method can be used for extracting first-order and total
  Sobol' indices, after the Shapley effects have been computed. It returns a list
  containing both sensitivity indices.
}

\references{
  Azadkia M., Chatterjee S. (2019). A simple measure of conditional dependence. 
  arXiv preprint arXiv:1910.12327.
  
		 Chatterjee, S., 2021, \emph{A new coefficient of correlation}, Journal of the American 
		 Statistical Association, 116:2009-2022. 
  
		 Gamboa, F., Gremaud, P., Klein, T., & Lagnoux, A., 2022, \emph{Global Sensitivity Analysis: 
		 a novel generation of mighty estimators based on rank statistics}, 
		 Bernoulli 28: 2345-2374.
    
  Broto B., Bachoc F. and Depecker M. (2020) \emph{Variance Reduction for Estimation
  of Shapley Effects and Adaptation to Unknown Input Distribution.} SIAM/ASA Journal
  on Uncertainty Quantification, 8(2).
  
  Castro J., Gomez D, Tejada J. (2009). \emph{Polynomial calculation of the Shapley value based 
  on sampling.} Computers & Operations Research, 36(5):1726-1730.
  
  M. Il Idrissi, V. Chabridon and B. Iooss (2021). \emph{Developments and applications
  of Shapley effects   to reliability-oriented sensitivity analysis with correlated inputs.}
  Environmental Modelling & Software, 143, 105115.
  
  M. Il Idrissi, V. Chabridon and B. Iooss (2021). \emph{Mesures d'importance relative  
  par decompositions de la performance de modeles de regression,} Preprint, 52emes Journees 
  de Statistiques de la Societe Francaise de Statistique (SFdS), pp. 497-502, 
  Nice, France, Juin 2021
}

\author{
Marouane Il Idrissi, Sebastien Da Veiga
}

\examples{
  \donttest{
  
library(parallel)
library(doParallel)
library(foreach)
library(gtools)
library(boot)
library(RANN)

###########################################################
# Linear Model with Gaussian correlated inputs

library(mvtnorm)

set.seed(1234)
n <- 1000
beta<-c(1,-1,0.5)
sigma<-matrix(c(1,0,0,
                0,1,-0.8,
                0,-0.8,1),
              nrow=3,
              ncol=3)

X <-rmvnorm(n, rep(0,3), sigma)
colnames(X)<-c("X1","X2", "X3")


y <- X\%*\%beta + rnorm(n,0,2)

# Without Bootstrap confidence intervals
x<-shapleysobol_knn(model=NULL, X=X,
            n.knn=3,
            noise=TRUE)
tell(x,y)
print(x)
plot(x)

#Using the extract method to get first-order and total Sobol' indices
extract(x)

# With Boostrap confidence intervals
x<-shapleysobol_knn(model=NULL, X=X,
            nboot=10, 
            n.knn=3,
            noise=TRUE,
            boot.level=0.7, 
            conf=0.95)
tell(x,y)
print(x)
plot(x)

#####################
# Extracting Sobol' indices with Bootstrap confidence intervals

nboot <- 10 # put nboot=50 for consistency

#Total Sobol' indices
x<-shapleysobol_knn(model=NULL, X=X,
            nboot=nboot, 
            n.knn=3,
            U=0,
            noise=TRUE,
            boot.level=0.7, 
            conf=0.95)
tell(x,y)
print(x)
plot(x)

#First-order Sobol' indices
x<-shapleysobol_knn(model=NULL, X=X,
            nboot=nboot, 
            n.knn=3,
            U=1,
            noise=TRUE,
            boot.level=0.7, 
            conf=0.95)
tell(x,y)
print(x)
plot(x)

#Closed Sobol' indices for specific subsets (list)
x<-shapleysobol_knn(model=NULL, X=X,
            nboot=nboot, 
            n.knn=3,
            U=list(c(1,2), c(1,2,3), 2),
            noise=TRUE,
            boot.level=0.7, 
            conf=0.95)
tell(x,y)
print(x)
plot(x)


#####################################################
# Test case: the non-monotonic Sobol g-function
# Example with a call to a numerical model
# First compute first-order indices with ranking
    
n <- 1000
X <- data.frame(matrix(runif(8 * n), nrow = n))
x <- shapleysobol_knn(model = sobol.fun, X = X, U = 1, method = "rank")
print(x)
plot(x)

library(ggplot2) ; ggplot(x)

# We can use the output sample generated for this estimation to compute total indices 
# without additional calls to the model
x2 <- shapleysobol_knn(model = NULL, X = X, U = 0, method = "knn", n.knn = 5)
tell(x2,x$y)
plot(x2)

ggplot(x2)


#####################################################
# Test case: the Ishigami function
# Example with given data and the use of approximate nearest neighbour search
n <- 5000
X <- data.frame(matrix(-pi+2*pi*runif(3 * n), nrow = n))
Y <- ishigami.fun(X)
x <- shapleysobol_knn(model = NULL, X = X, U = NULL, method = "knn", n.knn = 5, 
                       n.limit = 2000)
tell(x,Y)
plot(x)

library(ggplot2) ; ggplot(x)

# Extract first-order and total Sobol indices
x1 <- extract(x) ; print(x1)
    
######################################################
# Test case : Linear model (3 Gaussian inputs including 2 dependent) with scaling
# See Iooss and Prieur (2019)
library(mvtnorm) # Multivariate Gaussian variables
library(whitening) # For scaling
modlin <- function(X) apply(X,1,sum)
d <- 3
n <- 10000
mu <- rep(0,d)
sig <- c(1,1,2)
ro <- 0.9
Cormat <- matrix(c(1,0,0,0,1,ro,0,ro,1),d,d)
Covmat <- ( sig \%*\% t(sig) ) * Cormat
Xall <- function(n) mvtnorm::rmvnorm(n,mu,Covmat)
X <- Xall(n)
x <- shapleysobol_knn(model = modlin, X = X, U = NULL, method = "knn", n.knn = 5, 
                       rescale = TRUE, n.limit = 2000)
print(x)
plot(x)
}
}

\seealso{
\code{\link{sobolrank}}, \code{\link{sobolshap_knn}}, \code{\link{shapleyPermEx}}, \code{\link{shapleySubsetMc}}, \code{\link{lmg}}, \code{\link{pme_knn}}
}

\keyword{shapley}
