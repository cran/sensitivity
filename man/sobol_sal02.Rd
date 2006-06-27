\name{sobol.sal02}
\alias{sobol.sal02}
\alias{compute.sobol.sal02}
\alias{print.sobol.sal02}
\alias{plot.sobol.sal02}

\title{Sobol Non-linear Sensitivity Analysis}

\description{
  \code{sobol.sal02} is the implementation of the Sobol indices
  estimation method proposed by Saltelli. The particularity of this
  method is to estimate both first order and total order indices with
  the same sampling scheme. This function generates the design of
  experiments based on two independent samples and computes the
  requested indices.
}

\usage{
sobol.sal02(model = NULL, x1, x2,
            nboot = 0, conf = 0.95, ...)
\method{compute}{sobol.sal02}(sa, y = NULL)
}

\arguments{
  \item{model}{the model.}
  \item{x1}{the first random sample for factors.}
  \item{x2}{the second random sample for factors.}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level for bootstrap confidence intervals.}
  \item{sa}{the sensitivity analysis object.}
  \item{y}{the response.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
}

\details{
  \code{model} is a function or a predictor (a class with a
  \code{predict} method) computing the response \code{y} based on the
  sample given by \code{x}. If no model is specified, the indices will be
  computed when one gives the response.

  The Monte Carlo estimation requires two independent random samples
  \code{x1} and \code{x2}. They must have the same dimensions.

  The computing scheme is based on replacement of a column of a sample
  by the corresponding one of the other sample. We note
  \eqn{\{x1,i,x2\}}{{x1,i,x2}} the sample where all the columns are taken
  from x1, except the i-th one that is taken from x2.

  Then, the first order index corresponding to the i-th factor is given
  by:
  \deqn{S_i=\frac{D_i}{\mathrm{Var}(Y)}}{S_i = D_i / Var(Y)}
  where
  \deqn{
    D_i = \frac{1}{n-1} \sum_{j=1}^n f(x2_j) f(\{x1,i,x2\}_j) 
    - \frac{1}{n} \sum_{j=1}^n f(x1_j) f(x2_j)
  }{
    D_i = 1 / (n - 1) * sum( f(x2) * f({x1,i,x2}) ) 
    - 1 / n * sum( f(x1) * f(x2) )
  }

  and the total order index is given by:
  \deqn{St_i=1 - \frac{Dt_i}{\mathrm{Var}(Y)}}{St_i = 1 - Dt_i / Var(Y)}
  where
  \deqn{
    Dt_i = \frac{1}{n-1} \sum_{j=1}^n f(x1_j) f(\{x1,i,x2\}_j)
    - \left[ \frac{1}{n} \sum_{j=1}^n f(x1_j) \right]^2
  }{
    Dt_i = 1 / (n - 1) * sum( f(x1) * f({x1,i,x2}) )
    - [ 1 / n * sum( f(x1) )]^2
  }

  With this method, the estimation of all the first order indices and
  total indices require \eqn{n(p+2)}{n * (p + 2)} calls to the model,
  with \eqn{n}{n} the length of the samples \code{x1} and \code{x2},
  and \eqn{p}{p} the number of factors.
}

\value{
  \code{sobol.sal02} returns an object of class \code{"sobol.sal02"}.
  An object of class \code{"sobol.sal02"} is a list containing the following
  components:

  \item{model}{the model.}
  \item{x1}{the first random sample for factors.}
  \item{x2}{the second random sample for factors.}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level for bootstrap confidence intervals.}
  \item{x}{the factor sample.}
  \item{y}{the response used.}
  \item{S1}{the estimations of the first order indices.}
  \item{St}{the estimations of the total sensitivity indices.}
  \item{call}{the matched call.}
}

\references{
Saltelli, A., 2002, Making best use of model evaluations to compute
sensitivity indices. \emph{Computer Physics Communication}, \bold{145},
580--297.
}

\seealso{
  \code{\link{sensitivity}}
  \code{\link{compute}}
}

\examples{
# Test case : the non-monotonic Sobol g-function

# The method of sobol requires 2 samples
# There are 8 factors, all following the uniform distribution
# on [0,1]

n <- 3000
x1 <- data.frame(matrix(nr = n, nc = 8))
x2 <- data.frame(matrix(nr = n, nc = 8))
for (i in 1:8){
  x1[, i] <- runif(n)
  x2[, i] <- runif(n)
}

# sensitivity analysis

sa <- sobol.sal02(model = sobol.fun, x1 = x1, x2 = x2, nboot = 100)
print(sa)
plot(sa)
}

\keyword{misc}
