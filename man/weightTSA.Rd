\name{weightTSA}
\alias{weightTSA}

\title{
Weight-function to transform an output variable in order to perform Target Sensitivity Analysis (TSA)
}
\description{
Transformation function of one variable (vector sample)
}
\usage{
weightTSA(Y, c, upper = TRUE, type="indicTh", param=1)
}

\arguments{
  \item{Y}{The output vector}
  \item{c}{The threshold}
  \item{upper}{TRUE for upper threshold and FALSE for lower threshold}
  \item{type}{The weight function type ("indicTh", "zeroTh", logistic", "exp1side"):
  \itemize{
    \item indicTh : indicator-thresholding
    \item zeroTh : zero-thresholding (keeps the variable value above (upper=TRUE case) or below the threshold)
    \item logistic : logistic transformation at the threshold
    \item exp1side : exponential transformation above (upper=TRUE case) or below the threshold (see Raguet and Marrel)
    }}
  \item{param}{The parameter value for "logistic" and "exp1side" types}
}

\details{The weight functions depend on a threshold \eqn{c} and/or a smooth relaxation. These functions are defined as follows
  \deqn{if type = "indicTh" --> w = 1_(Y>c),}
  \deqn{if type = "zeroTh" --> w = Y*1_(Y>c),}
  \deqn{if type = "logistic" --> w = 1/(1 + exp(-param*(Y-c)/abs(c))),}
  \deqn{if type = "exp1side" --> w = exp{-max(c - Y, 0)/(param.\sigma(Y)/5)},}
  where \eqn{\sigma(Y)} is an estimation of the standard deviation of Y and \eqn{param = 1} is a parameter tuning the smoothness.
}

\value{
The vector sample of the transformed variable
}

\references{

H. Raguet and A. Marrel, \emph{Target and conditional sensitivity analysis with emphasis
on dependence measures}, Preprint, https://hal.archives-ouvertes.fr/hal-01694129

A. Spagnol, \emph{Kernel-based sensitivity indices for high-dimensional optimization problems},
PhD Thesis, Universite de Lyon, 2020

Spagnol A., Le Riche R., Da Veiga S. (2019), \emph{Global sensitivity analysis for optimization 
with variable selection}, SIAM/ASA J. Uncertainty Quantification, 7(2), 417--443.
}

\author{
B. Iooss
}

\examples{
n <- 100  # sample size
c <- 1.5
Y <- rnorm(n)
Yt <- weightTSA(Y, c)
}
