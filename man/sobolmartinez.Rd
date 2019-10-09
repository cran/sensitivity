\name{sobolmartinez}
\alias{sobolmartinez}
\alias{tell.sobolmartinez}
\alias{print.sobolmartinez}
\alias{plot.sobolmartinez}
\alias{ggplot.sobolmartinez}

\title{Monte Carlo Estimation of Sobol' Indices (formulas of Martinez (2011))}

\description{
  \code{sobolmartinez} implements the Monte Carlo estimation of
  the Sobol' indices for both first-order and total indices using 
  correlation coefficients-based formulas, at a total cost of 
  \eqn{(p+2) \times n}{(p + 2) * n} model evaluations. 
  These are called the Martinez estimators.
}

\usage{
sobolmartinez(model = NULL, X1, X2, nboot = 0, conf = 0.95, \dots)
\method{tell}{sobolmartinez}(x, y = NULL, return.var = NULL, \dots)
\method{print}{sobolmartinez}(x, \dots)
\method{plot}{sobolmartinez}(x, ylim = c(0, 1), y_col = NULL, y_dim3 = NULL, \dots)
\method{ggplot}{sobolmartinez}(x, ylim = c(0, 1), y_col = NULL, y_dim3 = NULL, \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{X1}{the first random sample.}
  \item{X2}{the second random sample.}
  \item{nboot}{the number of bootstrap replicates, or zero to use theoretical 
  formulas based on confidence interfaces of correlation coefficient 
  (Martinez, 2011).}
  \item{conf}{the confidence level for bootstrap confidence intervals.}
  \item{x}{a list of class \code{"sobol"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{return.var}{a vector of character strings giving further
    internal variables names to store in  the output object \code{x}.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{y_col}{an integer defining the index of the column of \code{x$y} to be
    used for plotting the corresponding sensitivity indices (only applies if 
    \code{x$y} is a matrix or an array). If set to \code{NULL} (as per default) 
    and \code{x$y} is a matrix or an array, the first column (respectively the 
    first element in the second dimension) of \code{x$y} is used (i.e. 
    \code{y_col = 1}).}
  \item{y_dim3}{an integer defining the index in the third dimension of 
    \code{x$y} to be used for plotting the corresponding sensitivity indices 
    (only applies if \code{x$y} is an array). If set to \code{NULL} (as per 
    default) and \code{x$y} is a three-dimensional array, the first element in 
    the third dimension of \code{x$y} is used (i.e. \code{y_dim3 = 1}).}
  \item{\dots}{for \code{sobolmartinez}: any other arguments for \code{model} 
    which are passed unchanged each time it is called}
}

\value{
  \code{sobolmartinez} returns a list of class \code{"sobolmartinez"}, 
  containing all the input arguments detailed before, plus the following 
  components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments.}
  \item{y}{either a vector, a matrix or a three-dimensional array of model 
    responses (depends on the output of \code{model}).}
  \item{V}{the estimations of normalized variances of the Conditional 
    Expectations (VCE) with respect to each factor and also with respect 
    to the complementary set of each factor ("all but \eqn{X_i}{Xi}").}
  \item{S}{the estimations of the Sobol' first-order indices.}
  \item{T}{the estimations of the Sobol' total sensitivity indices.}
  
  Users can ask more ouput variables with the argument
  \code{return.var} (for example, bootstrap outputs \code{V.boot},
  \code{S.boot} and \code{T.boot}).
}

\details{
This estimator supports missing values (NA or NaN) which can occur during the
simulation of the model on the design of experiments (due to code failure) 
even if Sobol' indices are no more rigorous variance-based sensitivity 
indices if missing values are present. In this case, a warning is displayed.

This version of \code{sobolmartinez} also supports matrices and 
three-dimensional arrays as output of \code{model}. Bootstrapping (including 
bootstrap confidence intervals) is also supported for matrix or array output. 
However, theoretical confidence intervals (for \code{nboot = 0}) are only 
supported for vector output. If the model output is a matrix or an array, 
\code{V}, \code{S} and \code{T} are matrices or arrays as well (depending on the
type of \code{y} and the value of \code{nboot}).

The bootstrap outputs \code{V.boot}, \code{S.boot} and \code{T.boot} can only be
returned if the model output is a vector (using argument \code{return.var}). For
matrix or array output, these objects can't be returned.
}

\references{
J-M. Martinez, 2011, \emph{Analyse de sensibilite globale par decomposition 
de la variance}, Presentation in the meeting of GdR Ondes and GdR MASCOT-NUM, 
January, 13th, 2011, Institut Henri Poincare, Paris, France.

M. Baudin, K. Boumhaout, T. Delage, B. Iooss and J-M. Martinez, 2016, Numerical stability of Sobol' indices estimation formula, Proceedings of the SAMO 2016 Conference, Reunion Island, France, December 2016
}

\author{
Bertrand Iooss, with contributions from Frank Weber (2016)
}

\seealso{
  \code{\link{sobol}, \link{sobol2002}, \link{sobolSalt}, \link{sobol2007}, \link{soboljansen}, \link{soboltouati}, \link{sobolMultOut}}
}

\examples{
# Test case : the non-monotonic Sobol g-function

# The method of sobol requires 2 samples
# There are 8 factors, all following the uniform distribution
# on [0,1]

library(boot)
n <- 1000
X1 <- data.frame(matrix(runif(8 * n), nrow = n))
X2 <- data.frame(matrix(runif(8 * n), nrow = n))

# sensitivity analysis

x <- sobolmartinez(model = sobol.fun, X1, X2, nboot = 0)
print(x)
plot(x)

library(ggplot2)
ggplot(x)

# Only for demonstration purposes: a model function returning a matrix
sobol.fun_matrix <- function(X){
  res_vector <- sobol.fun(X)
  cbind(res_vector, 2 * res_vector)
}
x_matrix <- sobolmartinez(model = sobol.fun_matrix, X1, X2)
plot(x_matrix, y_col = 2)
title(main = "y_col = 2")

# Also only for demonstration purposes: a model function returning a
# three-dimensional array
sobol.fun_array <- function(X){
  res_vector <- sobol.fun(X)
  res_matrix <- cbind(res_vector, 2 * res_vector)
  array(data = c(res_matrix, 5 * res_matrix), 
        dim = c(length(res_vector), 2, 2))
}
x_array <- sobolmartinez(model = sobol.fun_array, X1, X2)
plot(x_array, y_col = 2, y_dim3 = 2)
title(main = "y_col = 2, y_dim3 = 2")
}

\keyword{design}
