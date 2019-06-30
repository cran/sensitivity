\name{sobolTIIlo}
\alias{sobolTIIlo}
\alias{tell.sobolTIIlo}
\alias{print.sobolTIIlo}
\alias{plot.sobolTIIlo}
\alias{ggplot.sobolTIIlo}
\alias{plotFG.sobolTIIlo}

\title{Liu and Owen Estimation of Total Interaction Indices}

\description{\code{sobolTIIlo} implements the  asymptotically efficient formula of Liu and Owen (2006) for the estimation of total interaction indices as described e.g. in Section 3.4 of Fruth et al. (2014). Total interaction indices (TII) are superset indices of pairs of variables, thus give the total influence of each second-order interaction. The total cost of the method is \eqn{(1+N+\choose(N,2)) \times n} where \eqn{N}{N} is the number
   of indices to estimate. Asymptotic confidence intervals are provided. Via \code{plotFG} (which uses functions of the package \code{igraph}), the TIIs can be visualized in a so-called FANOVA graph as described in section 2.2 of Muehlenstaedt et al. (2012).}

\usage{
sobolTIIlo(model = NULL, X1, X2, conf = 0.95, \dots)
\method{tell}{sobolTIIlo}(x, y = NULL, \dots)
\method{print}{sobolTIIlo}(x, \dots)
\method{plot}{sobolTIIlo}(x, ylim = NULL, \dots)
\method{ggplot}{sobolTIIlo}(x, ylim = NULL, \dots)
\method{plotFG}{sobolTIIlo}(x)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{X1}{the first random sample.}
  \item{X2}{the second random sample.}
 \item{conf}{the confidence level for asymptotic confidence intervals, defaults to 0.95.}
  \item{x}{a list of class \code{"sobolTIIlo"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
  \item{ylim}{optional, the y limits of the plot.}
}

\value{
  \code{sobolTIIlo} returns a list of class \code{"sobolTIIlo"}, containing all
  the input arguments detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments.}
  \item{y}{a vector of model responses.}
  \item{V}{the estimation of the overall variance.}
  \item{tii.unscaled}{the unscaled estimations of the TIIs.}
  \item{tii.scaled}{the scaled estimations of the TIIs together with asymptotic confidence intervals.}
}

\references{
R. Liu, A. B. Owen, 2006,  \emph{Estimating mean dimensionality of analysis of variance decompositions}, JASA, 101 (474), 712--721.

J. Fruth, O. Roustant, S. Kuhnt, 2014, \emph{Total interaction index: A variance-based sensitivity index for second-order interaction screening}, J. Stat. Plan. Inference, 147, 212--223.

T. Muehlenstaedt, O. Roustant, L. Carraro, S. Kuhnt, 2012, \emph{Data-driven Kriging models based on FANOVA-decomposition}, Stat. Comput., 22 (3), 723--738.
}

\author{
Jana Fruth
}

\seealso{
  \code{\link{sobolTIIpf}}
}

\examples{
# Test case : the Ishigami function

# The method requires 2 samples
n <- 1000
X1 <- data.frame(matrix(runif(3 * n, -pi, pi), nrow = n))
X2 <- data.frame(matrix(runif(3 * n, -pi, pi), nrow = n))

# sensitivity analysis (the true values of the scaled TIIs are 0, 0.244, 0)
x <- sobolTIIlo(model = ishigami.fun, X1 = X1, X2 = X2)
print(x)

# plot of tiis and FANOVA graph
plot(x)

library(ggplot2)
ggplot(x)

\donttest{
library(igraph)
plotFG(x)
}


}

\keyword{design}
