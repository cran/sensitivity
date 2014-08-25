\name{sensiHSIC}
\alias{sensiHSIC}
\alias{tell.sensiHSIC}
\alias{print.sensiHSIC}
\alias{plot.sensiHSIC}

\title{Sensitivity Indices based on Hilbert-Schmidt Independence Criterion (HSIC)}

\description{ \code{sensiHSIC} conducts a sensitivity analysis where the impact of
an input variable is defined in terms of the distance between the input/output
joint probability distribution and the product of their marginals when they are
embedded in a Reproducing Kernel Hilbert Space (RKHS). This distance corresponds
to the Hilbert-Schmidt Independence Criterion (HSIC) proposed by Gretton et al.
(2005) and serves as a dependence measure between random variables, see Da Veiga
(2014) for an illustration in the context of sensitivity analysis.}

\usage{
  sensiHSIC(model = NULL, X, kernelX = "rbf", paramX = NA, 
            kernelY = "rbf", paramY = NA, nboot = 0, conf = 0.95, ...)
  \method{tell}{sensiHSIC}(x, y = NULL, \dots)
  \method{print}{sensiHSIC}(x, \dots)
  \method{plot}{sensiHSIC}(x, ylim = c(0, 1), \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
               defining the model to analyze.}
  \item{X}{a matrix or \code{data.frame} representing the input random sample.}
  \item{kernelX}{a string or a list of strings specifying the reproducing kernel
  to be used for the input variables. If only one kernel is provided, it is used
  for all input variables. Available choices are "rbf" (Gaussian), "laplace"
  (exponential), "dcov" (distance covariance, see details), "raquad" (rationale
  quadratic), "invmultiquad" (inverse multiquadratic), "linear" (Euclidean scalar
  product), "matern3" (Matern 3/2), "matern5" (Matern 5/2), "ssanova1" (kernel of
  Sobolev space of order 1) and "ssanova2" (kernel of Sobolev space of order 2).}
  \item{paramX}{a scalar or a vector of hyperparameters to be used in the input
  variable kernels. If only one scalar is provided, it is replicated for all input
  variables. By default \code{paramX} is equal to the standard deviation of the
  input variable for "rbf", "laplace", "raquad", "invmultiquad", "matern3" and
  "matern5" and to \code{1} for "dcov". Kernels "linear", "ssanova1" and "ssanova2"
  do not involve hyperparameters. If \code{kernelX} is a combination of kernels
  with and without hyperparameters, \code{paramX} must have a (dummy) value for the
  hyperparameter-free kernels, see examples below.}
  \item{kernelY}{a string specifying the reproducing kernel to be used for the
  output variable. Available choices are "rbf" (Gaussian), "laplace" (exponential),
  "dcov" (distance covariance, see details), "raquad" (rationale quadratic),
  "invmultiquad" (inverse multiquadratic), "linear" (Euclidean scalar product),
  "matern3" (Matern 3/2), "matern5" (Matern 5/2), "ssanova1" (kernel of Sobolev
  space of order 1) and "ssanova2" (kernel of Sobolev space of order 2).}
  \item{paramY}{a scalar to be used in the output variable kernel. By default
  \code{paramY} is equal to the standard deviation of the output variable for "rbf",
  "laplace", "raquad", "invmultiquad", "matern3" and "matern5" and to \code{1} for "dcov".
  Kernels "linear", "ssanova1" and "ssanova2" do not involve hyperparameters.}
  \item{nboot}{the number of bootstrap replicates}
  \item{conf}{the confidence level for confidence intervals.}
  \item{x}{a list of class \code{"sensiHSIC"} storing the state of the
           sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{\dots}{any other arguments for \code{model} which are passed
               unchanged each time it is called.}
}

\value{
  \code{sensiHSIC} returns a list of class \code{"sensiHSIC"}, containing all
  the input arguments detailed before, plus the following components:
    
    \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments.}
  \item{y}{a vector of model responses.}
  \item{S}{the estimations of HSIC sensitivity indices.}
}

\details{
  The HSIC sensitivity indices are obtained as a normalized version of the Hilbert-Schmidt independence criterion: \deqn{S_i^{HSIC} = \frac{HSIC(X_i,Y)}{\sqrt{HSIC(X_i,X_i)}\sqrt{HSIC(Y,Y)}},}{%
  Si = HSIC(Xi,Y) / (\sqrt HSIC(Xi,Xi) \sqrt HSIC(Y,Y)),} see Da Veiga (2014) for details.
  When \code{kernelX="dcov"} and \code{kernelY="dcov"}, the kernel is given by \eqn{k(u,u')=1/2(||u||+||u'||-||u-u'||)} and the sensitivity index is equal to the distance correlation introduced by Szekely et al. (2007) as was recently proven by Sejdinovic et al. (2013).
}
\references{
  
  Da Veiga S. (2014), \emph{Global sensitivity analysis with dependence measures},
  Journal of Statistical Computation and Simulation, in press. 
  \url{http://hal.archives-ouvertes.fr/hal-00903283}
  
  Gretton A., Bousquet O., Smola A., Scholkopf B. (2005), 
  \emph{Measuring statistical dependence with hilbert-schmidt norms}, 
  Jain S, Simon H, Tomita E, editors: Algorithmic learning theory, 
  Lecture Notes in Computer Science, Vol. 3734, Berlin: Springer, 63--77.
  
  Sejdinovic D., Sriperumbudur B., Gretton A., Fukumizu K., (2013), 
  \emph{Equivalence of distance-based and RKHS-based statistics in hypothesis
  testing}, Annals of Statistics 41(5), 2263--2291.
  
  Szekely G.J., Rizzo M.L., Bakirov N.K. (2007), 
  \emph{Measuring and testing dependence by correlation of distances}, 
  Annals of Statistics 35(6), 2769--2794.
}

\author{
  Sebastien Da Veiga, Snecma
}

\seealso{
\code{\link[ks]{kde}, \link{sensiFdiv}}
}

\examples{
 \dontrun{ 
  # Test case : the non-monotonic Sobol g-function
  # Only one kernel is provided with default hyperparameter value
  n <- 100
  X <- data.frame(matrix(runif(8 * n), nrow = n))
  x <- sensiHSIC(model = sobol.fun, X, kernelX = "raquad", kernelY = "rbf")
  print(x)
  
  # Test case : the Ishigami function
  # A list of kernels is given with default hyperparameter value
  n <- 100
  X <- data.frame(matrix(-pi+2*pi*runif(3 * n), nrow = n))
  x <- sensiHSIC(model = ishigami.fun, X, kernelX = c("rbf","matern3","dcov"), 
                  kernelY = "rbf")
  print(x)
  
  # A combination of kernels is given and a dummy value is passed for 
  # the first hyperparameter
  x <- sensiHSIC(model = ishigami.fun, X, kernelX = c("ssanova1","matern3","dcov"), 
                  paramX = c(1,2,1), kernelY = "ssanova1")
  print(x)
 }
}
