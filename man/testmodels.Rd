\name{testmodels}
\alias{testmodels}
\alias{sobol.fun}
\alias{ishigami.fun}
\alias{morris.fun}

\title{Test Models For Sensitivity Analysis}

\description{
  These functions are standard testcase for sensitivity
  analysis benchmarks. There are: the g-function of Sobol, the function
  of Ishigami and the function of Morris.
}

\usage{
sobol.fun(x)
ishigami.fun(x)
morris.fun(x)
}

\arguments{
  \item{x}{the \code{matrix} or \code{data.frame} containing the input
    values.}
}

\value{
  All these functions return a \code{numeric} vector containig the values
  of the function.
}

\references{
Saltelli, A., Chan, K. and Scott, E. M., 2000, \emph{Sensitivity
  analysis}. Wiley, 33--45.
}

\seealso{
  \code{\link{sensitivity}}
}

\examples{
# Scatterplots for Sobol g-function versus X1 and X2

n <- 2000
x <- matrix(nr = n, nc = 8)
for (i in 1:8)
  x[, i] <- runif(n)

y <- sobol.fun(x)

plot(x[, 1], y)
plot(x[, 2], y)

# Scatterplots for the Ishigami function versus X1, X2 and X3

n <- 2000
x <- matrix(nr = n, nc = 3)
for (i in 1:3)
  x[, i] <- runif(n, min = -pi, max = pi)

y <- ishigami.fun(x)

plot(x[, 1], y)
plot(x[, 2], y)
plot(x[, 3], y)
}

\keyword{misc}
