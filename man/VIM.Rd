\name{VIM}
\alias{VIM}
\alias{print.VIM}
\alias{plot.VIM}
\alias{ggplot.VIM}

\title{Variance-based Importance Measures in linear model}

\description{
  \code{VIM} summarizes some linear variance-based importance measures useful 
  in data analysis/machine learning context (dependent inputs' case): 
  VIF (i.e. variance inflation factor which is a multicollinearity metric), squared SRC,
  squared PCC, LMG and PMVD, as well as the R2 and Q2 of the linear regression model
}

\usage{
VIM(X, y, logistic = FALSE, nboot = 0, 
    conf = 0.95, max.iter = 1000, parl = NULL)
\method{print}{VIM}(x, \dots)
\method{plot}{VIM}(x, ylim = c(0,1), \dots)
\method{ggplot}{VIM}(data, mapping = aes(), \dots, ylim = c(0,1), 
  environment = parent.frame())
}

\arguments{
  \item{X}{a matrix or data frame containing the observed covariates
  (i.e., features, input variables...).}
  \item{y}{a numeric vector containing the observed outcomes (i.e.,
  dependent variable). If \code{logistic=TRUE}, can be a numeric vector
  of zeros and ones, or a logical vector, or a factor.}
  \item{logistic}{logical. If \code{TRUE}, the analysis is done via a
  logistic regression(binomial GLM).}
  \item{nboot}{the number of bootstrap replicates for the computation
  of confidence intervals.}
  \item{conf}{the confidence level of the bootstrap confidence intervals.}
  \item{max.iter}{if \code{logistic=TRUE}, the maximum number of iterative 
  optimization steps allowed for the logistic regression. Default is \code{1000}.} 
  \item{parl}{number of cores on which to parallelize the computation. If
  \code{NULL}, then no parallelization is done.}
  \item{x}{the object returned by \code{VIM}.}
  \item{data}{the object returned by \code{VIM}.}
  \item{ylim}{the y-coordinate limits of the plot.}
  \item{mapping}{Default list of aesthetic mappings to use for plot. If not specified, 
    must be supplied in each layer added to the plot.}
  \item{environment}{[Deprecated] Used prior to tidy evaluation.}
  \item{\dots}{arguments to be passed to methods, such as graphical
    parameters (see \code{par}).}
}

\value{
  \code{VIM} returns a list of class \code{"VIM"}, containing the following
  components:

  \item{call}{the matched call.}
  \item{R2}{a data frame containing the estimations of the R2.}
  \item{Q2}{a data frame containing the estimations of the Q2.}
  \item{VIF}{a data frame containing the estimations of the VIF.}
  \item{SRC2}{a data frame containing the estimations of the squared SRC.}
  \item{PCC2}{a data frame containing the estimations of the squared PCC.}
  \item{LMG}{a data frame containing the estimations of the LMG.}
  \item{PMVD}{a data frame containing the estimations of the PMVD.}
  \item{X}{the observed covariates.}
  \item{y}{the observed outcomes.}
  \item{logistic}{logical. \code{TRUE} if the analysis has been made by
  logistic regression.}
  \item{nboot}{number of bootstrap replicates.}
  \item{max.iter}{if \code{logistic=TRUE}, the maximum number of iterative 
  optimization steps allowed for the logistic regression. Default is \code{1000}.} 
  \item{parl}{number of chosen cores for the computation.}
  \item{conf}{level for the confidence intervals by bootstrap.}
}

\details{
  This function cannot be used with categorical inputs.

  For logistic regression (\code{logistic=TRUE}), the \eqn{R^2}{R-squared}
  value is equal to:
  \deqn{R^2 = 1-\frac{\textrm{model deviance}}{\textrm{null deviance}}}{R-squared
  = 1 - (model deviance)/(null deviance)}

  If too many cores for the machine are passed on to the \code{parl} argument,
  the chosen number of cores is defaulted to the available cores minus one.
}

\references{
  L. Clouvel, B. Iooss, V. Chabridon, M. Il Idrissi and F. Robin, 2024,
  \emph{An overview of variance-based importance measures in the linear regression context: 
  comparative analyses and numerical tests}, Socio-Environmental Systems Modelling, vol. 7, 
  18681, 2025, doi:10.18174/sesmo.1868.
  \url{https://hal.science/hal-04102053}
}

\author{
Bertrand Iooss
}

\examples{
\donttest{

library(parallel)
library(boot)
library(car)

library(mvtnorm)

set.seed(1234)
n <- 100
sigma<-matrix(c(1,0,0,0.9, 0,1,-0.8,0, 0,-0.8,1,0, 0.9,0,0,1), nr=4, nc=4)

############################
# Gaussian correlated inputs

X <- as.data.frame(rmvnorm(n, rep(0,4), sigma))
colnames(X) <- c("X1","X2","X3","X4")

#############################
# Linear Model with small noise, two correlated inputs (X2 and X3) and 
# one dummy input (X4) correlated with another (X1)
epsilon <- rnorm(n,0,0.1)
y <- with(X, X1 - X2 + 0.5 * X3 + epsilon)

# Without Bootstrap confidence intervals
x <- VIM(X, y)
print(x)
plot(x)
library(ggplot2) ; ggplot(x)

# With Boostrap confidence intervals
x <- VIM(X, y, nboot=100, conf=0.9)
print(x)
plot(x)
library(ggplot2) ; ggplot(x)

############################
# Logistic Regression (same regression model)

epsilon <- rnorm(n,0,0.1)
y <- with(X, X1 - X2 + 0.5 * X3 + epsilon > 0)

x <- VIM(X, y, logistic = TRUE)
print(x)
plot(x)
library(ggplot2) ; ggplot(x)
}
}

\seealso{
\code{\link{src}}, \code{\link{pcc}}, \code{\link{src}}, \code{\link{lmg}}, \code{\link{pmvd}}
}

\keyword{regression}
