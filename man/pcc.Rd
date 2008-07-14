\name{pcc}
\alias{pcc}
\alias{print.pcc}
\alias{plot.pcc}

\title{Partial Correlation Coefficients}

\description{
  \code{pcc} computes the Partial Correlation Coefficients (PCC), or
  Partial Rank Correlation Coefficients (PRCC), which are sensitivity
  indices based on linear (resp. monotonic) assumptions, in the case of
  (linearly) correlated factors.
}

\usage{
pcc(X, y, rank = FALSE, nboot = 0, conf = 0.95)
\method{print}{pcc}(x, \dots)
\method{plot}{pcc}(x, ylim = c(-1,1), \dots)
}

\arguments{
  \item{X}{a data frame (or object coercible by \code{as.data.frame})
    containing the design of experiments (model input variables).}
  \item{y}{a vector containing the responses corresponding to the design
    of experiments (model output variables).}
  \item{rank}{logical. If \code{TRUE}, the analysis is done on the
    ranks.}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level of the bootstrap confidence intervals.}
  \item{x}{the object returned by \code{pcc}.}
  \item{ylim}{the y-coordinate limits of the plot.}
  \item{\dots}{arguments to be passed to methods, such as graphical
    parameters (see \code{par}).}
}

\details{
}

\value{
  \code{pcc} returns a list of class \code{"pcc"}, containing the following
  components:

  \item{call}{the matched call.}
  \item{PCC}{a data frame containing the estimations of the PCC
    indices, bias and confidence intervals (if \code{rank = TRUE}).}
  \item{PRCC}{a data frame containing the estimations of the PRCC
    indices, bias and confidence intervals (if \code{rank = TRUE}).}
}

\references{
A. Saltelli, K. Chan and E. M. Scott eds, 2000, \emph{Sensitivity
  Analysis}, Wiley.
}

\examples{
# a 100-sample with X1 ~ U(0.5, 1.5)
#                   X2 ~ U(1.5, 4.5)
#                   X3 ~ U(4.5, 13.5)
n <- 100
X <- data.frame(X1 = runif(n, 0.5, 1.5),
                X2 = runif(n, 1.5, 4.5),
                X3 = runif(n, 4.5, 13.5))

# linear model : Y = X1 + X2 + X3
y <- with(X, X1 + X2 + X3)

# sensitivity analysis
x <- pcc(X, y, nboot = 100)
print(x)
#plot(x) # TODO: find another example...
}

\seealso{
\code{\link{src}}
}

\keyword{regression}
