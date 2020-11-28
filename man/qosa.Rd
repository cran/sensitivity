\name{qosa}
\alias{qosa}
\alias{tell.qosa}
\alias{print.qosa}
\alias{plot.qosa}
\alias{ggplot.qosa}

\title{Quantile-oriented sensitivity analysis}

\description{ \code{qosa} implements the estimation of first-order quantile-oriented sensitivity indices 
as defined in Fort et al. (2016) with a kernel-based estimator of conditonal probability density functions 
closely related to the one proposed by Maume-Deschamps and Niang (2018). 
\code{qosa} also supports a kernel-based estimation of Sobol first-order indices (i.e. Nadaraya-Watson).
}

\usage{
qosa(model = NULL, X1, X2 = NULL, type = "quantile", alpha = 0.1, split.sample = 2/3, 
nsample = 1e4, nboot = 0, conf = 0.95, \dots)
\method{tell}{qosa}(x, y = NULL, \dots)
\method{print}{qosa}(x, \dots)
\method{plot}{qosa}(x, ylim = c(0, 1), \dots)
\method{ggplot}{qosa}(x, ylim = c(0, 1), \dots)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method, defining the model to analyze.}
  \item{X1}{a random sample of the inputs used for the estimation of conditional probability density functions. 
  If \code{X2} is NULL, \code{X1} is split in two samples, with the first \code{split.sample} proportion of observations assigned to \code{X1} and the rest to \code{X2}.}
  \item{X2}{a random sample of the inputs used to evaluate the conditional probability density functions. 
  If NULL, it is constructed with the last \code{(1-split.sample)} proportion of observations from \code{X1}, see above.}
  \item{type}{a string specifying which first-order sensitivity indices must be estimated: quantile-oriented indices (\code{type="quantile"}) 
  or Sobol' indices (\code{type="mean"}).}
  \item{alpha}{if \code{type="quantile"} the quantile level.}
  \item{split.sample}{if \code{X2=NULL} the proportion of observations from \code{X1} assigned 
  to the estimation of conditional probability density functions.}
  \item{nsample}{the number of samples from the conditional probability density functions used 
  to estimate the conditional quantiles (if \code{type="quantile"}) or the conditional means (if \code{type="mean"}).}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level for confidence intervals.}
  \item{x}{a list of class \code{"sobolrank"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
}

\value{
  \code{qosa} returns a list of class \code{"qosa"}, containing all
  the input arguments detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments.}
  \item{X1}{a \code{data.frame} containing the design of experiments used for the estimation of conditional probability density functions.}
  \item{X}{a \code{data.frame} containing the design of experiments used for the evaluation of conditional probability density functions.}
  \item{y}{a vector of model responses.}
  \item{S}{the estimations of the Sobol' sensitivity indices.}
}
\details{
		 Quantile-oriented sensitivty indices were defined as a special case of sensitivity indices based on contrast functions in Fort et al. (2016).
		 The estimator used by \code{qosa} follows closely the one proposed by Maume-Deschamps & Niang (2018). 
		 The only difference is that Maume-Deschamps and Niang (2018) use the following kernel-based estimate of the conditional cumulative distribution function:
		 \deqn{\hat{F}(y\Vert X=x) = \frac{ \sum_{i=1}^n K_{h_x}(x - X_i) \bold{1}\{Y_i< y\}}{\sum_{i=1}^n K_{h_x}(x - X_i)}}
		 whereas we use
		 \deqn{\hat{F}(y\vert X=x) = \frac{ \sum_{i=1}^n K_{h_x}(x - X_i) \int_{-\infty}^y K_{h_y}(t - Y_i)dt} {\sum_{i=1}^n K_{h_x}(x - X_i)},}
		 meaning that \eqn{\bold{1}\{Y_i< y\}} is replaced by \eqn{\int_{-\infty}^y K_{h_y}(t - Y_i)dt = \Phi (\frac{y-Y_i}{h_y})} where \eqn{\Phi} 
		 is the cumulative distribution function of the standard normal distribution (since kernel \eqn{K} is Gaussian). 
		 The two definitions thus coincide when \eqn{h_y \rightarrow 0}. Our formula arises from a kernel density estimator of the joint pdf with a diagonal bandwidth. 
		 In a future version, it will be genralized to a general bandwidth matrix for improved performance.

}
\references{
		 Fort, J. C., Klein, T., and Rachdi, N. (2016). New sensitivity analysis subordinated to a contrast. Communications in Statistics-Theory and Methods, 45(15), 4349-4364.
		 
		 Maume-Deschamps, V., and Niang, I. (2018). Estimation of quantile oriented sensitivity indices. Statistics & Probability Letters, 134, 122-127.
}
\author{
Sebastien Da Veiga
}
\examples{
 \donttest{
library(ks)
library(ggplot2)
library(boot)

# Test case : difference of two exponential distributions (Fort et al. (2016))
# We use two samples with different size
n1 <- 5000
X1 <- data.frame(matrix(rexp(2 * n1,1), nrow = n1))
n2 <- 1000
X2 <- data.frame(matrix(rexp(2 * n2,1), nrow = n2))
Y1 <- X1[,1] - X1[,2]
Y2 <- X2[,1] - X2[,2]
x <- qosa(model = NULL, X1, X2, type = "quantile", alpha = 0.1)
tell(x,c(Y1,Y2))
print(x)
ggplot(x)

# Test case : difference of two exponential distributions (Fort et al. (2016))
# We use only one sample
n <- 10000
X <- data.frame(matrix(rexp(2 * n,1), nrow = n))
Y <- X[,1] - X[,2]
x <- qosa(model = NULL, X1 = X, type = "quantile", alpha = 0.7)
tell(x,Y)
print(x)
ggplot(x)

# Test case : the Ishigami function
# We estimate first-order Sobol' indices (by specifying 'mean')
n <- 5000
X <- data.frame(matrix(-pi+2*pi*runif(3 * n), nrow = n))
x <- qosa(model = ishigami.fun, X1 = X, type = "mean", nboot = 50)
print(x)
ggplot(x)
}
}