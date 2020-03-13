\name{sensiHSIC}
\alias{sensiHSIC}
\alias{tell.sensiHSIC}
\alias{print.sensiHSIC}
\alias{plot.sensiHSIC}
\alias{ggplot.sensiHSIC}

\title{Sensitivity Indices based on Hilbert-Schmidt Independence Criterion (HSIC)}

\description{ \code{sensiHSIC} conducts a sensitivity analysis where the impact of
an input variable is defined in terms of the distance between the input/output
joint probability distribution and the product of their marginals when they are
embedded in a Reproducing Kernel Hilbert Space (RKHS). This distance corresponds
to the Hilbert-Schmidt Independence Criterion (HSIC) proposed by Gretton et al.
(2005) and serves as a dependence measure between random variables, see Da Veiga
(2015) for an illustration in the context of sensitivity analysis. The use of universal
kernels ensures equivalence between HSIC nullity and independence of X and Y. In this
case, a statistical test of independence with HSIC measure as statistic can be built.
H0 is "X and Y are independent", against H1: X and Y are dependent. P-value can be
computed either with asymptotic approximation (Gamma approximation), either with
permutation method. See Meynaoui et al. (2019) for details. \code{sensiHSIC} can also 
be used to perform a Target Sensitivity Analysis (TSA). The basic idea of TSA is to 
measure the influence of the inputs on a critical domain of the model output, see 
Spagnol et al. (2019), Raguet and Marrel (2018) for details.
}

\usage{
  sensiHSIC(model = NULL, X, target = NULL, kernelX = "rbf", paramX = NA,
            kernelY = "rbf", paramY = NA, nboot = 0, conf = 0.95,
            estimator.type = "V-stat", test.method = "Asymptotic", B = 5000, 
            crit.option = list(stop.criterion = "screening", alpha = 0.05, Bstart = 100,
            Bfinal = 5000, Bbatch = 100, lo = 100), ...)
  \method{tell}{sensiHSIC}(x, y = NULL, \dots)
  \method{print}{sensiHSIC}(x, \dots)
  \method{plot}{sensiHSIC}(x, ylim = c(0, 1), \dots)
  \method{ggplot}{sensiHSIC}(x, ylim = c(0, 1), \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
               defining the model to analyze.}
  \item{X}{a matrix or \code{data.frame} representing the input random sample.}
  \item{target}{list of parameters to perform Target Sensitivity Analysis (TSA).
  \code{c}: threshold. \code{upper}: TRUE for upper threshold and FALSE for lower threshold.
  \code{type}: the weight function type ("indicTh", "logistic", "exp1side").
  \code{param}: the parameter value for "logistic" and "exp1side" types.}
  \item{kernelX}{a string or a list of strings specifying the reproducing kernel
  to be used for the input variables. If only one kernel is provided, it is used
  for all input variables. Available choices are "rbf" (Gaussian), "laplace"
  (exponential), "dcov" (distance covariance, see details), "raquad" (rationale
  quadratic), "invmultiquad" (inverse multiquadratic), "linear" (Euclidean scalar
  product), "matern3" (Matern 3/2), "matern5" (Matern 5/2), "ssanova1" (kernel of
  Sobolev space of order 1) and "ssanova2" (kernel of Sobolev space of order 2)}
  \item{paramX}{a scalar or a vector of hyperparameters to be used in the input
  variable kernels. If only one scalar is provided, it is replicated for all input
  variables. By default \code{paramX} is equal to the standard deviation of the
  input variable for "rbf", "laplace", "raquad", "invmultiquad", "matern3" and
  "matern5" and to \code{1} for "dcov". Kernels "linear", "ssanova1" and "ssanova2"
  do not involve hyperparameters. If \code{kernelX} is a combination of kernels
  with and without hyperparameters, \code{paramX} must have a (dummy) value for the
  hyperparameter-free kernels, see examples below}
  \item{kernelY}{a string specifying the reproducing kernel to be used for the
  output variable. Available choices are "rbf" (Gaussian), "laplace" (exponential),
  "dcov" (distance covariance, see details), "raquad" (rationale quadratic),
  "invmultiquad" (inverse multiquadratic), "linear" (Euclidean scalar product),
  "matern3" (Matern 3/2), "matern5" (Matern 5/2), "ssanova1" (kernel of Sobolev
  space of order 1) and "ssanova2" (kernel of Sobolev space of order 2).}
  \item{paramY}{a scalar to be used in the output variable kernel. By default
  \code{paramY} is equal to the standard deviation of the output variable for "rbf",
  "laplace", "raquad", "invmultiquad", "matern3" and "matern5" and to \code{1} for "dcov".
  Kernels "linear", "ssanova1" and "ssanova2" do not involve hyperparameters}
  \item{nboot}{the number of bootstrap replicates}
  \item{conf}{the confidence level for confidence intervals}
  \item{estimator.type}{a string specifying the type of HSIC estimator. Two types of
  estimators are available, V-statistic of U-statistic. If estimator.type = "V-stat"
  (default value), the HSIC is estimated with a biased (but asymptotically unbiased)
  estimator, more practical for numerical implementation. If estimator.type = "U-stat",
  the unbiased estimator is used. The variance is of order o(1/n) for both estimators
  (n being the sample size, i.e. number of rows of X). More details in Meynaoui et al., 2019}
  \item{test.method}{a string specifying the method used to compute the p-value of the
  HSIC-based independence test. If test.method = "Asymptotic" (default value), an asympotic
  approximation with Gamma distribution is used. If test.method = "Permutation", a
  permutation method based on B boostrap samples is used to estimate the p-value in a
  non-asymptotic framework. Permutation-test are recommended for low sample size n.
  More details in Meynaoui et al., 2019. If test.method = "Seq_Permutation", an iterative
  permutation method is used to estimate the p-value. This approach bypasses the a-priori 
  choice of the   number of permutations (B) and the sequential estimation is stopped 
  according   to the user's objective (see crit.option for details).}
  \item{B}{number of boostrap samples used in the permutation method for the estimation of
  P-values in independence test. B is used only if test.method is "Permutation"}
  \item{crit.option}{a list of parameters used only if test.method = "Seq_Permutation".
  \code{stop.criterion}:"ranking" or "screening". \code{alpha}: significance level of the test.
  \code{Bstart}: initial number of boostrap samples used for P-values estimation. \code{Bfinal}:
  final number   of boostrap samples. \code{Bbatch}: batch of boostrap samples used in the 
  iterative estimation.
  \code{lo}: parameter depending on the stability we want to achieve.}
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
  \item{X}{a \code{data.frame} containing the design of experiments}
  \item{y}{a vector of model responses}
  \item{S}{the estimations of normalized HSIC sensitivity indices (also denoted R2HSIC)}
  \item{HSICXY}{the estimation of HSIC sensitivity indices (numerator in S formula)}
  \item{Pvalue}{the estimation of P-value of the independence test based on HSIC statistic}

}

\details{
  The HSIC sensitivity indices are obtained as a normalized version of the Hilbert-Schmidt
  independence criterion:
  \deqn{S_i^{HSIC} = R^2_{HSIC} = \frac{HSIC(X_i,Y)}{\sqrt{HSIC(X_i,X_i)}\sqrt{HSIC(Y,Y)}},}{%
  Si = HSIC(Xi,Y) / (\sqrt HSIC(Xi,Xi) \sqrt HSIC(Y,Y)),} see Da Veiga (2014) for details.
  When \code{kernelX="dcov"} and \code{kernelY="dcov"}, the kernel is given by
  \eqn{k(u,u')=1/2(||u||+||u'||-||u-u'||)} and the sensitivity index is equal to the distance
  correlation introduced by Szekely et al. (2007) as was recently proven by Sejdinovic et al. (2013).
  The target sensitivity measures are defined via weight functions \eqn{w} which depend on 
  \eqn{c = threshold}. 
  Indicator function \eqn{1_c} and smooth relaxations are available (according to target$type):
  \deqn{if type = "indicTh" --> w = 1_(Y>c),}
  \deqn{if type = "logistic" --> w = 1/(1 + exp(-param*(Y-c)/abs(c))), cf. Spagnol et al. (2019)}
  \deqn{if type = "exp1side" --> w = exp{-max(c - Y, 0)/(param.\sigma(Y)/5)}, 
  cf. Raguet and Marrel (2018)}   where \eqn{\sigma(Y)} is an estimation of the standard deviation 
  of Y and \eqn{param = 1} is a parameter tuning the smoothness.
}

\references{

  Da Veiga S. (2015), \emph{Global sensitivity analysis with dependence measures},
  Journal of Statistical Computation and Simulation, 85(7), 1283--1305.

  Gretton A., Bousquet O., Smola A., Scholkopf B. (2005),
  \emph{Measuring statistical dependence with hilbert-schmidt norms},
  Jain S, Simon H, Tomita E, editors: Algorithmic learning theory,
  Lecture Notes in Computer Science, Vol. 3734, Berlin: Springer, 63--77.

  Meynaoui A., Marrel A., and Laurent B. (2019).
  \emph{New statistical methodology for second level global sensitivity analysis},
  Preprint, ArXiv 1902.07030.

  Raguet H., Marrel A. (2018).
  \emph{Target and conditional sensitivity analysis with emphasis on dependence measures},
  Preprint, ArXiv 1801.10047.

  Sejdinovic D., Sriperumbudur B., Gretton A., Fukumizu K., (2013),
  \emph{Equivalence of distance-based and RKHS-based statistics in hypothesis
  testing}, Annals of Statistics 41(5), 2263--2291.

  Spagnol A., Le Riche R., Da Veiga S. (2019),
  \emph{Global sensitivity analysis for optimization with variable selection},
  SIAM/ASA J. Uncertainty Quantification, 7(2), 417--443.

  Szekely G.J., Rizzo M.L., Bakirov N.K. (2007),
  \emph{Measuring and testing dependence by correlation of distances},
  Annals of Statistics 35(6), 2769--2794.
}

\author{
  Sebastien Da Veiga, Amandine Marrel, Anouar Meynaoui, Reda El Amri
}

\seealso{
\code{\link[ks]{kde}, \link{sensiFdiv}}
}

\examples{
 \donttest{
  library(ggplot2)
  library(boot)

  # Test case : the non-monotonic Sobol g-function
  # Only one kernel is provided with default hyperparameter value

  n <- 100
  X <- data.frame(matrix(runif(8 * n), nrow = n))

  # HSIC-based GSA and asymptotic independence test
  x <- sensiHSIC(model = sobol.fun, X, kernelX = "raquad", kernelY = "rbf",
    nboot = 100, test.method = "Asymptotic")
  print(x)

  # HSIC-based GSA and permutation independence test
  x <- sensiHSIC(model = sobol.fun, X, kernelX = "rbf", kernelY = "rbf",
    estimator.type = "U-stat", test.method = "Permutation")
  print(x)

  # HSIC-based GSA and independence test with optimized number of permutations
  x <- sensiHSIC(model = sobol.fun, X, kernelX = "raquad", kernelY = "rbf",
    nboot = 100, test.method = "Seq_Permutation",
    crit.option = list(stop.criterion = "ranking", alpha = 0.05, Bstart = 100, Bfinal = 3000, 
    Bbatch = 100, lo = 100))
  print(x)

  # HSIC-Target GSA in case of given model
  x <- sensiHSIC(model = sobol.fun, X, target = list(c = 0.4, upper = TRUE, type = "indicTh", 
    param = 1), kernelX = "raquad", kernelY = "rbf",nboot = 100, test.method = "Permutation")

  # HSIC-target GSA in case of given data
  x <- sensiHSIC(model = NULL, X, target = list(c = 0.4, upper = TRUE, type = "indicTh", param = 1),
    kernelX = "raquad", kernelY = "rbf",nboot = 100, test.method = "Permutation")
  y <- sobol.fun(X)
  tell(x,y)

  # Test case : the Ishigami function

  # A list of kernels is given with default hyperparameter value
  n <- 100
  X <- data.frame(matrix(-pi+2*pi*runif(3 * n), nrow = n))
  x <- sensiHSIC(model = ishigami.fun, X, kernelX = c("rbf","matern3","dcov"),
                  kernelY = "rbf")
  print(x)
  ggplot(x)

  # A combination of kernels is given and a dummy value is passed for
  # the first hyperparameter
  x <- sensiHSIC(model = ishigami.fun, X, kernelX = c("ssanova1","matern3","dcov"),
                  paramX = c(1,2,1), kernelY = "ssanova1")
  print(x)
  ggplot(x)

  # Example in case of given data
  n <- 100
  X <- data.frame(matrix(-pi+2*pi*runif(3 * n), nrow = n))
  Y <- ishigami.fun(X)

  x <- sensiHSIC(model = NULL, X)
  tell(x,Y)
  print(x)
  ggplot(x)

 }
}
