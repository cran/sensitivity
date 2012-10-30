\name{sb}
\alias{sb}
\alias{ask.sb}
\alias{tell.sb}
\alias{print.sb}
\alias{plot.sb}

\title{Sequential Bifurcations}

\description{ \code{sb} implements the Sequential Bifurcations screening
  method (Bettonvil and Kleijnen 1996). \strong{This is an alpha version
    that might strongly evolve in the future}.}

\usage{
sb(p, sign = rep("+", p), interaction = FALSE)
\method{ask}{sb}(x, i = NULL, \dots)
\method{tell}{sb}(x, y, \dots)
\method{print}{sb}(x, \dots)
\method{plot}{sb}(x, \dots)
}

\arguments{
  \item{p}{number of factors.}
  \item{sign}{a vector fo length \code{p} filled with \code{"+"} and
  \code{"-"}, giving the (assumed) signs of the factors effects.}
  \item{interaction}{a boolean, \code{TRUE} if the model is supposed to
    be with interactions, \code{FALSE} otherwise.}
  \item{x}{a list of class \code{"sb"} storing the state of the
    screening study at the current iteration.}
  \item{y}{a vector of model responses.}
  \item{i}{an integer, used to force a wanted bifurcation instead of that
    proposed by the algorithm.}
  \item{\dots}{not used.}
}

\details{
  The model without interaction is
  \deqn{Y=\beta_0 + \sum_{i=1}^p \beta_i X_i}{
    Y = beta_0 + sum_{i=1}^p beta_i X_i}
  while the model with interactions is
  \deqn{Y=\beta_0 + \sum_{i=1}^p \beta_i X_i + \sum_{1 \leq i < j \leq
    p} \gamma_{ij} X_i X_j}{
    Y = beta_0 + sum_{i=1}^p beta_i X_i + sum_{1 <= i < j <= p} gamma_{ij} X_i X_j}
  In both cases, the factors are assumed to be uniformly distributed on
  \eqn{[-1,1]}{[-1,1]}. This is a difference with Bettonvil
  et al. where the factors vary across \eqn{[0,1]}{[0,1]} in the former
  case, while \eqn{[-1,1]}{[-1,1]} in the latter.

  Another difference with Bettonvil et al. is that in the current
  implementation, the groups are splitted right in the middle.
}

\value{
  \code{sb} returns a list of class \code{"sb"}, containing all
  the input arguments detailed before, plus the following components:
  
  \item{i}{the vector of bifurcations.}
  \item{y}{the vector of observations.}
  \item{ym}{the vector of mirror observations (model with interactions
    only).}
  
  The groups effects can be displayed with the \code{print} method.
}

\references{
B. Bettonvil and J. P. C. Kleijnen, 1996, \emph{Searching for important
factors in simulation models with many factors: sequential
bifurcations}, European Journal of Operational Research, 96, 180--194.
}

\examples{
# a model with interactions
p <- 50
beta <- numeric(length = p)
beta[1:5] <- runif(n = 5, min = 10, max = 50)
beta[6:p] <- runif(n = p - 5, min = 0, max = 0.3)
beta <- sample(beta)
gamma <- matrix(data = runif(n = p^2, min = 0, max = 0.1), nrow = p, ncol = p)
gamma[lower.tri(gamma, diag = TRUE)] <- 0
gamma[1,2] <- 5
gamma[5,9] <- 12
f <- function(x) { return(sum(x * beta) + (x \%*\% gamma \%*\% x))}

# 10 iterations of SB
sa <- sb(p, interaction = TRUE)
for (i in 1 : 10) {
  x <- ask(sa)
  y <- list()
  for (i in names(x)) {
    y[[i]] <- f(x[[i]])
  }
  tell(sa, y)
}
print(sa)
plot(sa)
}

\keyword{design}
