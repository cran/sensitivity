\name{sobolTIIpf}
\alias{sobolTIIpf}
\alias{tell.sobolTIIpf}
\alias{print.sobolTIIpf}
\alias{plot.sobolTIIpf}
\alias{ggplot.sobolTIIpf}
\alias{plotFG.sobolTIIpf}
\alias{plotFG}

\title{Pick-freeze Estimation of Total Interaction Indices}

\description{\code{sobolTIIpf} implements the pick-freeze estimation of total interaction indices as described in Section 3.3 of Fruth et al. (2014). Total interaction indices (TII) are superset indices of pairs of variables, thus give the total influence of each second-order interaction. The pick-freeze estimation enables the strategy to reuse evaluations of Saltelli (2002). The total costs are \eqn{(1+N) \times n} where \eqn{N}{N} is the number of indices to estimate. Via \code{plotFG}, the TIIs can be visualized in a so-called FANOVA graph as described in section 2.2 of Muehlenstaedt et al. (2012).}

\usage{
sobolTIIpf(model = NULL, X1, X2, \dots)
\method{tell}{sobolTIIpf}(x, y = NULL, \dots)
\method{print}{sobolTIIpf}(x, \dots)
\method{plot}{sobolTIIpf}(x, ylim = NULL, \dots)
\method{ggplot}{sobolTIIpf}(data, mapping = aes(), ylim = NULL, \dots, environment
                 = parent.frame())
\method{plotFG}{sobolTIIpf}(x)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{X1}{the first random sample.}
  \item{X2}{the second random sample.}
  \item{x}{a list of class \code{"sobolTIIpf"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{data}{a list of class \code{"sobolTIIpf"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{mapping}{Default list of aesthetic mappings to use for plot. If not specified, 
    must be supplied in each layer added to the plot.}
  \item{environment}{[Deprecated] Used prior to tidy evaluation.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
  \item{ylim}{optional, the y limits of the plot.}
}

\value{
  \code{sobolTIIpf} returns a list of class \code{"sobolTIIpf"}, containing all
  the input arguments detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments.}
  \item{y}{a vector of model responses.}
  \item{V}{the estimation of the overall variance.}
  \item{tii.unscaled}{the unscaled estimations of the TIIs together.}
  \item{tii.scaled}{the scaled estimations of the TIIs.}
}

\references{
J. Fruth, O. Roustant, S. Kuhnt, 2014, \emph{Total interaction index: A variance-based sensitivity index for second-order interaction screening}, J. Stat. Plan. Inference, 147, 212--223.

A. Saltelli, 2002, \emph{Making best use of model evaluations to compute sensitivity indices}, Comput. Phys. Commun., 145, 580-297.

T. Muehlenstaedt, O. Roustant, L. Carraro, S. Kuhnt, 2012, \emph{Data-driven Kriging models based on FANOVA-decomposition}, Stat. Comput., 22 (3), 723--738.
}

\author{
Jana Fruth
}

\seealso{
  \code{\link{sobolTIIlo}}
}

\examples{
# Test case : the Ishigami function

# The method requires 2 samples
n <- 1000
X1 <- data.frame(matrix(runif(3 * n, -pi, pi), nrow = n))
X2 <- data.frame(matrix(runif(3 * n, -pi, pi), nrow = n))

# sensitivity analysis (the true values are 0, 0.244, 0)
x <- sobolTIIpf(model = ishigami.fun, X1 = X1, X2 = X2)
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
