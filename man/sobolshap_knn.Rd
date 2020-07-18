\name{sobolshap_knn}
\alias{sobolshap_knn}
\alias{tell.sobolshap_knn}
\alias{extract.sobolshap_knn}
\alias{print.sobolshap_knn}
\alias{plot.sobolshap_knn}
\alias{ggplot.sobolshap_knn}

\title{Flexible sensitivity analysis via ranking / nearest neighbours}

\description{\code{sobolshap_knn} implements the estimation of several sensitivity indices using 
only N model evaluations via ranking (following Gamboa et al. (2020) and Chatterjee (2019)) 
or nearest neighbour search (Broto et al. (2020) and Azadkia & Chatterjee (2020)). 
It can be used with categorical inputs (which are transformed with one-hot encoding), 
dependent inputs and multiple outputs. Sensitivity indices of any group of inputs can be computed,
which means that in particular first-order/total Sobol indices and Shapley effects are accessible. 
For large sample sizes, the nearest neightbour algorithm can be significantly accelerated 
by using approximate nearest neighbour search.
}

\usage{
  sobolshap_knn(model = NULL, X, id.cat = NULL, U = NULL, method = "knn", n.knn = 2, 
                return.shap = FALSE, rescale = FALSE, n.limit = 2000, noise = FALSE, \dots)
  \method{tell}{sobolshap_knn}(x, y = NULL, \dots)
  \method{extract}{sobolshap_knn}(x, \dots)
  \method{print}{sobolshap_knn}(x, \dots)
  \method{plot}{sobolshap_knn}(x, ylim = c(0, 1), type.multout = "lines", \dots)
  \method{ggplot}{sobolshap_knn}(x, ylim = c(0, 1), type.multout = "lines", \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method, defining the model to analyze.}
  \item{X}{a random sample of the inputs.}
  \item{id.cat}{a vector with the indices of the categorical inputs.}
  \item{U}{an integer equal to 0 (total Sobol indices) or 1 (first-order Sobol indices) 
  or a list of vector indices defining the subsets of inputs whose sensitivity indices 
  must be computed or a matrix of 0s and 1s where each row encodes a subset of inputs 
  whose sensitivity indices must be computed (see examples) or NULL (all possible subsets).}
  \item{method}{the algorithm to be used for estimation, either "rank" or "knn", see details.}
  \item{n.knn}{the number of nearest neighbours used for estimation if \code{method="rank"}.}
  \item{return.shap}{a logical indicating if Shapley effects must be estimated, 
  can only be TRUE if \code{U=NULL}.}
  \item{rescale}{a logical indicating if continuous inputs must be rescaled before distance computations.
  If TRUE, continuous inputs are first whitened with the ZCA-cor whitening procedure 
  (cf. whiten() function in package \code{whitening}). If the inputs are independent, 
  this first step will have a very limited impact. Then, the resulting whitened inputs 
  are individually modified via a copula transform such that each input has the same scale.}
  \item{n.limit}{the sample size limit above which approximate nearest neighbour search is activated, 
  only used if \code{method="rank"}.}
  \item{noise}{a logical which is TRUE if the model or the output sample is noisy, see details.}
  \item{x}{a list of class \code{"sobolshap_knn"} storing the state of the sensitivity study 
  (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{type.multout}{the plotting method in the case of multiple outputs, either "points" or "lines", 
  see examples.}
  \item{\dots}{any other arguments for \code{model} which are passed unchanged each time it is called.}
}

\value{
  \code{sobolshap_knn} returns a list of class \code{"sobolrank"}, containing all
  the input arguments detailed before, plus the following components:
    
    \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments.}
  \item{y}{a vector of model responses.}
  \item{U}{the subsets of inputs for which sensitivity indices have been computed.}
  \item{S}{the estimations of the Sobol sensitivity indices (see details).}
  \item{Shap}{the estimations of Shapley effects, if return.shap was set to TRUE.}
  \item{order}{0 (total indices), 1 (first-order indices) or NULL. Used for plotting defaults.}
}

\details{
    For \code{method="rank"}, the estimator is defined in Gamboa et al. (2020) 
    following Chatterjee (2019). For first-order indices it is based on an input ranking 
    (same algorithm as in \code{sobolrank}) while for higher orders, 
    it uses an approximate heuristic solution of the traveling salesman problem 
    applied to the input sample distances (cf. TSP() function in package \code{TSP}).
    For \code{method="knn"}, ranking and TSP are replaced by a nearest neighbour search 
    as proposed in Broto et al. (2020) and in Azadkia & Chatterjee (2020) for a similar coefficient. 
    The algorithm is the same as in \code{shapleySubsetMc} but with an optimized implementation. 
    In particular, the distance used for subsets with mixed inputs (continuous and categorical) 
    are the same but here the additional one-hot encoding of categorical variables makes it possible to
    work only with Euclidean distances. Furthermore, a fast approximate nearest neighbour search is also
    available, which is strongly recommended for large sample sizes. The main difference 
    with \code{shapleySubsetMc} is that here we use the entire N sample to compute all indices, 
    while in \code{shapleySubsetMc} the user can specify a total cost Ntot which performs 
    a specific allocation of sample sizes to the estimation of each index. 
    In addition, the \code{weights} option is not available here yet.
    If the outputs are noisy, the argument \code{noise} can be used: it only has an impact on the
    estimation of one specific sensitivity index, namely \eqn{Var(E(Y|X1,\ldots,Xp))/Var(Y)}. 
    If there is no noise this index is equal to 1, while in the presence of noise it must be estimated.
    The \code{extract} method is useful if in a first step the Shapley effects have been computed 
    and thus sensitivity indices for all possible subsets are available. 
    The resulting \code{sobolshap_knn} object can be post-treated by \code{extract} 
    to get first-order and total Sobol indices very easily.
}

\references{
    Azadkia M., Chatterjee S. (2019). A simple measure of conditional dependence. 
    arXiv preprint arXiv:1910.12327.
    Broto B., Bachoc F., Depecker M. (2020), Variance reduction for estimation of Shapley effects 
    and adaptation to unknown input distribution, SIAM/ASA Journal of Uncertainty Quantification, 
    8:693-716.
    Chatterjee S. (2019). A new coefficient of correlation. arXiv preprint arXiv:1909.10140.
    Gamboa F., Gremaud P., Klein T., Lagnoux A. (2020). Global Sensitivity Analysis: 
    a new generation of mighty estimators based on rank statistics. arXiv preprint arXiv:2003.01772.
}

\author{
Sebastien Da Veiga
}

\seealso{
    \code{\link{sobolrank}, \link{shapleySubsetMc}}
}

\examples{
  \donttest{
    # Test case: the non-monotonic Sobol g-function
    # Example with a call to a numerical model
    # First compute first-order indices with ranking
    n <- 1000
    X <- data.frame(matrix(runif(8 * n), nrow = n))
    x <- sobolshap_knn(model = sobol.fun, X = X, U = 1, method = "rank")
    print(x)
    library(ggplot2)
    ggplot(x)
    # We can use the output sample generated for this estimation to compute total indices 
    # without additional calls to the model
    x2 <- sobolshap_knn(model = NULL, X = X, U = 0, method = "knn", n.knn = 5)
    tell(x2,x$y)
    ggplot(x2)
    
    # Test case: the Ishigami function
    # Example with given data and the use of approximate nearest neighbour search
    library(RANN)
    n <- 5000
    X <- data.frame(matrix(-pi+2*pi*runif(3 * n), nrow = n))
    Y <- ishigami.fun(X)
    x <- sobolshap_knn(model = NULL, X = X, U = NULL, method = "knn", n.knn = 5, 
                       return.shap = TRUE, n.limit = 2000)
    tell(x,Y)
    library(ggplot2)
    ggplot(x)
    # We can also extract first-order and total Sobol indices
    x1 <- extract(x)
    print(x1)
    
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
    x <- sobolshap_knn(model = modlin, X = X, U = NULL, method = "knn", n.knn = 5, 
                       return.shap = TRUE, rescale = TRUE, n.limit = 2000)
    print(x)
    
    # Test case: functional toy fct 'Arctangent temporal function'
    n <- 3000
    X <- data.frame(matrix(runif(2*n,-7,7), nrow = n))
    Y <- atantemp.fun(X)
    x <- sobolshap_knn(model = NULL, X = X, U = NULL, method = "knn", n.knn = 5, 
                       return.shap = TRUE, n.limit = 2000)
    tell(x,Y)
    library(ggplot2)
    library(reshape2)
    ggplot(x, type.multout="lines")
    }
}
