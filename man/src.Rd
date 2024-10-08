\name{src}
\alias{src}
\alias{print.src}
\alias{plot.src}
\alias{ggplot.src}

\title{Standardized Regression Coefficients}

\description{
  \code{src} computes the Standardized Regression Coefficients
  (SRC), or the Standardized Rank Regression Coefficients (SRRC), which
  are sensitivity indices based on linear or monotonic assumptions in
  the case of independent factors.
}

\usage{
src(X, y, rank = FALSE, logistic = FALSE, nboot = 0, conf = 0.95)
\method{print}{src}(x, \dots)
\method{plot}{src}(x, ylim = c(-1,1), \dots)
\method{ggplot}{src}(data,  mapping = aes(), ylim = c(-1, 1), \dots, environment
                 = parent.frame())
}

\arguments{
  \item{X}{a data frame (or object coercible by \code{as.data.frame})
    containing the design of experiments (model input variables).}
  \item{y}{a vector containing the responses corresponding to the design
    of experiments (model output variables).}
  \item{rank}{logical. If \code{TRUE}, the analysis is done on the
    ranks.}
  \item{logistic}{logical. If \code{TRUE}, the analysis is done via a
    logistic regression (binomial GLM).}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level of the bootstrap confidence intervals.}
  \item{x}{the object returned by \code{src}.}
  \item{data}{the object returned by \code{src}.}
  \item{ylim}{the y-coordinate limits of the plot.}
  \item{mapping}{Default list of aesthetic mappings to use for plot. If not specified, 
    must be supplied in each layer added to the plot.}
  \item{environment}{[Deprecated] Used prior to tidy evaluation.}
  \item{\dots}{arguments to be passed to methods, such as graphical
    parameters (see \code{par}).}
}

\value{
  \code{src} returns a list of class \code{"src"}, containing the following
  components:

  \item{call}{the matched call.}
  \item{SRC}{a data frame containing the estimations of the SRC
    indices, bias and confidence intervals (if \code{rank = FALSE}).}
  \item{SRRC}{a data frame containing the estimations of the SRRC
    indices, bias and confidence intervals (if \code{rank = TRUE}).}
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
  
A. Saltelli, K. Chan and E. M. Scott eds, 2000, \emph{Sensitivity
  Analysis}, Wiley.
  
}

\author{
Gilles Pujol and Bertrand Iooss
}

\examples{

# a 100-sample with X1 ~ U(0.5, 1.5)
#                   X2 ~ U(1.5, 4.5)
#                   X3 ~ U(4.5, 13.5)

library(boot)
n <- 100
X <- data.frame(X1 = runif(n, 0.5, 1.5),
                X2 = runif(n, 1.5, 4.5),
                X3 = runif(n, 4.5, 13.5))

# linear model : Y = X1 + X2 + X3

y <- with(X, X1 + X2 + X3)

# sensitivity analysis

x <- src(X, y, nboot = 100)
print(x)
plot(x)

library(ggplot2)
ggplot(x)
}

\seealso{
\code{\link{pcc}}, \code{\link{lmg}}, \code{\link{pmvd}}
}

\keyword{regression}
