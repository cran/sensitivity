\name{pme_knn}
\alias{pme_knn}
\alias{tell.pme_knn}
\alias{print.pme_knn}
\alias{plot.pme_knn}
\alias{ggplot.pme_knn}

\title{Data-given proportional marginal effects estimation via nearest-neighbors procedure}

\description{
\code{pme_knn} computes the proportional marginal effects (PME), from Herin et al. (2024) 
via a nearest neighbor estimation.
Parallelized computations are possible to accelerate the estimation process.
  It can be used with categorical inputs (which are transformed with one-hot encoding before 
  computing the nearest-neighbors), dependent inputs and multiple outputs.
  For large sample sizes, the nearest neighbour algorithm can be significantly accelerated 
  by using approximate nearest neighbour search.
}

\usage{
pme_knn(model=NULL, X, method = "knn", tol = NULL, marg = T, n.knn = 2, 
          n.limit = 2000, noise = F, rescale = F, nboot = NULL, 
          boot.level = 0.8, conf=0.95, parl=NULL, \dots)
\method{tell}{pme_knn}(x, y, \dots)
\method{print}{pme_knn}(x, \dots)
\method{plot}{pme_knn}(x, ylim = c(0,1), \dots)
\method{ggplot}{pme_knn}(data, mapping = aes(), ylim = c(0, 1), \dots, environment
                 = parent.frame())
}

\arguments{
  \item{model}{a function defining the model to analyze, taking X as an argument.}
  \item{X}{a matrix or data frame containing the observed inputs.}
  \item{method}{the algorithm to be used for estimation, either "rank" or "knn",
  see details. Default is \code{method="knn"}.}
  \item{tol}{tolerance under which an input is considered as being a zero input. 
  See details.}
  \item{marg}{whether to chose the closed Sobol' (\code{FALSE}) or total Sobol' 
  (\code{TRUE}) indices as value functions.}
  \item{n.knn}{the number of nearest neighbours used for estimation.}
  \item{n.limit}{sample size limit above which approximate nearest neighbour
  search is activated.}
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
  \item{x}{the object returned by \code{pme_knn}.}
  \item{data}{the object returned by \code{pme_knn}.}
  \item{y}{a numeric univariate vector containing the observed outputs.}
  \item{ylim}{the y-coordinate limits for plotting.}
  \item{mapping}{Default list of aesthetic mappings to use for plot. If not specified, 
    must be supplied in each layer added to the plot.}
  \item{environment}{[Deprecated] Used prior to tidy evaluation.}
  \item{\dots}{additional arguments to be passed to \code{model}, or to the 
  methods, such as graphical parameters (see \code{par}).}
}

\value{
  \code{pme_knn} returns a list of class \code{"pme_knn"}:

  \item{call}{the matched call.}
  \item{PME}{the estimations of the PME indices.}
  \item{VE}{the estimations of the closed Sobol' indices for all possible sub-models.}
  \item{indices}{list of all subsets corresponding to the structure of VE.}
  \item{method}{which estimation method has been used.}
  \item{conf_int}{a matrix containing the estimations, biais and confidence
  intervals by bootstrap (if \code{nboot>0}).}
  \item{X}{the observed covariates.}
  \item{y}{the observed outcomes.}
  \item{n.knn}{value of the \code{n.knn} argument.}
  \item{rescale}{wheter the design matrix has been rescaled.}
  \item{n.limit}{value of the \code{n.limit} argument.}
  \item{boot.level}{value of the \code{boot.level} argument.}
  \item{noise}{wheter the PME must sum up to one or not.}
  \item{boot}{logical, wheter bootstrap confidence interval estimates have 
  been performed.}
  \item{nboot}{value of the \code{nboot} argument.}
  \item{parl}{value of the \code{parl} argument.}
  \item{conf}{value of the \code{conf} argument.}
  \item{marg}{value of the \code{marg} argument.}
  \item{tol}{value of the \code{tol} argument.}
  
}

\details{
  For \code{method="rank"}, the estimator is defined in Gamboa et al. (2020) 
  following Chatterjee (2019).For first-order indices it is based on an input 
  ranking   (same algorithm as in \code{sobolrank}) while for higher orders, 
  it uses an approximate heuristic solution of the traveling salesman problem 
  applied to the input sample distances (cf. TSP() function in package 
  \code{TSP}).  For \code{method="knn"}, ranking and TSP are replaced by a 
  nearest neighbour search   as proposed in Broto et al. (2020) and in Azadkia 
  & Chatterjee (2020) for a similar coefficient. 

  The computation is done using the subset procedure, defined in Broto, Bachoc
  and Depecker (2020), that is computing all the Sobol' closed indices for all 
  possible sub-models first, and then computing the proportional values 
  recursively, as detailed in Feldman (2005), but using an extension to 
  non strictly positive games (Herin et al., 2024).

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
  
  If \code{marg = TRUE} (default), the chosen value function to compute the 
  proportional values are the total Sobol' indices (dual of the underlying 
  cooperative game). If \code{marg = FALSE}, then the closed Sobol' indices 
  are used instead. Differences may appear between the two.
  
  Zero inputs are defined by the \code{tol} argument. If \code{null},
  then inputs with:
  \deqn{S^T_{\{i\}}) = 0}{Total Sobol' index=0}
  are considered as zero input in the detection of spurious variables. If 
  provided, zero inputs are detected when:
  \deqn{S^T_{\{i\}} \leq \textrm{tol}}{Total Sobol' index <= tol}
}

\references{
  Azadkia M., Chatterjee S., 2021), \emph{A simple measure of conditional dependence}, 
  Ann. Statist. 49(6):3070-3102.
  
		 Chatterjee, S., 2021, \emph{A new coefficient of correlation}, Journal of the American 
		 Statistical Association, 116:2009-2022. 
  
		 Gamboa, F., Gremaud, P., Klein, T., & Lagnoux, A., 2022, \emph{Global Sensitivity Analysis: 
		 a novel generation of mighty estimators based on rank statistics}, 
		 Bernoulli 28: 2345-2374.
    
  Broto B., Bachoc F. and Depecker M. (2020) \emph{Variance Reduction for Estimation
  of Shapley Effects and Adaptation to Unknown Input Distribution.} SIAM/ASA Journal
  on Uncertainty Quantification, 8(2).
  
  M. Herin, M. Il Idrissi, V. Chabridon and B. Iooss, Proportional marginal effects for 
  sensitivity analysis with correlated inputs, Proceedings of the 10th International
  Conferenceon Sensitivity Analysis of Model Output (SAMO 2022), p 42-43, 
  Tallahassee, Florida, March 2022.
  
  M. Herin, M. Il Idrissi, V. Chabridon and B. Iooss, \emph{Proportional marginal effects 
  for global sensitivity analysis}, SIAM/ASA Journal of Uncertainty Quantification, 
  12:667-692 2024
  
  M. Il Idrissi, V. Chabridon and B. Iooss (2021). \emph{Developments and applications
  of Shapley effects   to reliability-oriented sensitivity analysis with correlated inputs.}
  Environmental Modelling & Software, 143, 105115.
  
  B. Iooss, V. Chabridon and V. Thouvenot, \emph{Variance-based importance 
  measures for machine learning model interpretability}, Congres lambda-mu23,
  Saclay, France, 10-13 octobre 2022
  \url{https://hal.science/hal-03741384}
  
  Feldman, B. (2005) \emph{Relative Importance and Value} SSRN Electronic Journal.
  
}

\author{
Marouane Il Idrissi, Margot Herin
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
x<-pme_knn(model=NULL, X=X,
            n.knn=3,
            noise=TRUE)
tell(x,y)
print(x)
plot(x)

# With Boostrap confidence intervals
x<-pme_knn(model=NULL, X=X,
            nboot=10, 
            n.knn=3,
            noise=TRUE,
            boot.level=0.7, 
            conf=0.95)
tell(x,y)
print(x)
plot(x)

#####################################################
# Test case: the Ishigami function
# Example with given data and the use of approximate nearest neighbour search
n <- 5000
X <- data.frame(matrix(-pi+2*pi*runif(3 * n), nrow = n))
Y <- ishigami.fun(X)
x <- pme_knn(model = NULL, X = X,  method = "knn", n.knn = 5, 
                       n.limit = 2000)
tell(x,Y)
plot(x)

library(ggplot2) ; ggplot(x)

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
x <- pme_knn(model = modlin, X = X, method = "knn", n.knn = 5, 
                       rescale = TRUE, n.limit = 2000)
print(x)
plot(x)
}
}

\seealso{
\code{\link{sobolrank}}, \code{\link{shapleysobol_knn}}, \code{\link{shapleyPermEx}}, \code{\link{shapleySubsetMc}}, \code{\link{lmg}}, \code{\link{pmvd}}
}

\keyword{proportional values}
