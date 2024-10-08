\name{pcc}
\alias{pcc}
\alias{print.pcc}
\alias{plot.pcc}
\alias{ggplot.pcc}

\title{Partial Correlation Coefficients}

\description{
  \code{pcc} computes the Partial Correlation Coefficients (PCC),
  Semi-Partial Correlation Coefficients (SPCC), Partial Rank Correlation 
  Coefficients (PRCC) or Semi-Partial Rank Correlation Coefficients (SPRCC), 
  which are variance-based measures based on linear (resp. monotonic) 
  assumptions, in the case of (linearly) correlated factors.
}

\usage{
pcc(X, y, rank = FALSE, semi = FALSE, logistic = FALSE, nboot = 0, conf = 0.95)
\method{print}{pcc}(x, \dots)
\method{plot}{pcc}(x, ylim = c(-1,1), \dots)
\method{ggplot}{pcc}(data, mapping = aes(), \dots, environment
                 = parent.frame(), ylim = c(-1,1))
}

\arguments{
  \item{X}{a data frame (or object coercible by \code{as.data.frame})
    containing the design of experiments (model input variables).}
  \item{y}{a vector containing the responses corresponding to the design
    of experiments (model output variables).}
  \item{rank}{logical. If \code{TRUE}, the analysis is done on the
    ranks.}
  \item{semi}{logical. If \code{TRUE}, semi-PCC are computed.}
  \item{logistic}{logical. If \code{TRUE}, the analysis is done via a
    logistic regression (binomial GLM).}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level of the bootstrap confidence intervals.}
  \item{x}{the object returned by \code{pcc}.}
  \item{data}{the object returned by \code{pcc}.}
  \item{ylim}{the y-coordinate limits of the plot.}
  \item{mapping}{Default list of aesthetic mappings to use for plot. If not specified, 
    must be supplied in each layer added to the plot.}
  \item{environment}{[Deprecated] Used prior to tidy evaluation.}
  \item{\dots}{arguments to be passed to methods, such as graphical
    parameters (see \code{par}).}
}

\value{
  \code{pcc} returns a list of class \code{"pcc"}, containing the following
  components:

  \item{call}{the matched call.}
  \item{PCC}{a data frame containing the estimations of the PCC
    indices, bias and confidence intervals (if \code{rank = TRUE} 
    and \code{semi = FALSE}).}
  \item{PRCC}{a data frame containing the estimations of the PRCC
    indices, bias and confidence intervals (if \code{rank = TRUE}
    and \code{semi = FALSE}).}
  \item{SPCC}{a data frame containing the estimations of the PCC
    indices, bias and confidence intervals (if \code{rank = TRUE} 
    and \code{semi = TRUE}).}
  \item{SPRCC}{a data frame containing the estimations of the PRCC
    indices, bias and confidence intervals (if \code{rank = TRUE}
    and \code{semi = TRUE}).}
}

\details{
  Logistic regression model (\code{logistic = TRUE}) and rank-based indices
  (\code{rank = TRUE}) are incompatible.
}

\references{
  
  L. Clouvel, B. Iooss, V. Chabridon, M. Il Idrissi and F. Robin, 2023,
  \emph{An overview of variance-based importance measures in the linear regression context: 
  comparative analyses and numerical tests}, Preprint.
  \url{https://hal.science/hal-04102053}

B. Iooss, V. Chabridon and V. Thouvenot, \emph{Variance-based importance 
measures for machine learning model interpretability}, Congres lambda-mu23,
Saclay, France, 10-13 octobre 2022
\url{https://hal.science/hal-03741384}
  
J.W. Johnson and J.M. LeBreton, 2004, \emph{History and use of relative 
importance indices in organizational research}, Organizational 
Research Methods, 7:238-257.

A. Saltelli, K. Chan and E. M. Scott eds, 2000, \emph{Sensitivity
  Analysis}, Wiley.
}

\author{
Gilles Pujol and Bertrand Iooss
}

\examples{
\donttest{
# a 100-sample with X1 ~ U(0.5, 1.5)
#                   X2 ~ U(1.5, 4.5)
#                   X3 ~ U(4.5, 13.5)
library(boot)
n <- 100
X <- data.frame(X1 = runif(n, 0.5, 1.5),
                X2 = runif(n, 1.5, 4.5),
                X3 = runif(n, 4.5, 13.5))

# linear model : Y = X1^2 + X2 + X3
y <- with(X, X1^2 + X2 + X3)

# sensitivity analysis
x <- pcc(X, y, nboot = 100)
print(x)
plot(x)

library(ggplot2)
ggplot(x)
ggplot(x, ylim = c(-1.5,1.5))

x <- pcc(X, y, semi = TRUE, nboot = 100)
print(x)
plot(x)
}
}

\seealso{
\code{\link{src}}, \code{\link{lmg}}, \code{\link{pmvd}}
}

\keyword{regression}
