\name{src}
\alias{src}
\alias{print.src}
\alias{plot.src}

\title{Standardized Regression Coefficients}

\description{
  \code{src} computes the Standardized Regression Coefficients
  (SRC), or the Standardized Rank Regression Coefficients (SRRC), which
  are sensitivity indices based on linear or monotonic assumptions in
  the case of independent factors.
}

\usage{
src(X, y, rank = FALSE, nboot = 0, conf = 0.95)
\method{print}{src}(x, \dots)
\method{plot}{src}(x, ylim = c(-1,1), \dots)
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
  \item{x}{the object returned by \code{src}.}
  \item{ylim}{the y-coordinate limits of the plot.}
  \item{\dots}{arguments to be passed to methods, such as graphical
    parameters (see \code{par}).}
}

\details{
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

x <- src(X, y, nboot = 100)
print(x)
plot(x)
}

\seealso{
\code{\link{pcc}}
}

\keyword{regression}
