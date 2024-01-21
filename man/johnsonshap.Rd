\name{johnsonshap}
\alias{johnsonshap}
\alias{print.johnsonshap}
\alias{plot.johnsonshap}
\alias{ggplot.johnsonshap}

\title{Johnson-Shapley indices}

\description{
  \code{johnsonshap} computes the Johnson-Shapley indices for correlated input relative importance. 
  These indices allocate a share of the output variance to each input based on the relative 
  weight allocation (RWA) system, in the case of dependent or correlated inputs.
  WARNING: This function does not work yet.
}

\usage{
johnsonshap(model = NULL, X, N, rank = FALSE, nboot = 0, conf = 0.95)
\method{print}{johnsonshap}(x, \dots)
\method{plot}{johnsonshap}(x, ylim = c(0,1), \dots)
\method{ggplot}{johnsonshap}(data,  mapping = aes(), ylim = c(0, 1), \dots, environment
                 = parent.frame())
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{X}{a data frame (or object coercible by \code{as.data.frame})
    containing a design of experiments (model input variables).}
  \item{N}{an integer giving the size of each replicated design for
  the Sobol' indices computations via the sobolrep() fct.}
  \item{rank}{logical. If \code{TRUE}, the analysis is done on the
    ranks.}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level of the bootstrap confidence intervals.}
  \item{x}{the object returned by \code{johnsonshap}.}
  \item{data}{the object returned by \code{johnsonshap}.}
  \item{ylim}{the y-coordinate limits of the plot.}
  \item{mapping}{Default list of aesthetic mappings to use for plot. If not specified, 
    must be supplied in each layer added to the plot.}
  \item{environment}{[Deprecated] Used prior to tidy evaluation.}
  \item{\dots}{arguments to be passed to methods, such as graphical
    parameters (see \code{par}).}
}

\value{
  \code{johnsonshap} returns a list of class \code{"johnsonshap"}, containing the following
  components:

  \item{call}{the matched call.}
  \item{johshap}{a data frame containing the estimations of the Johnson-Shapley
    indices, bias and confidence intervals.}
}

\details{
  No
}

\references{

B. Iooss and L. Clouvel, \emph{Une methode d'approximation des effets 
de Shapley en grande dimension}, 54emes  Journees de Statistique, 
Bruxelles, Belgique, July 3-7, 2023

}

\author{
Bertrand Iooss
}

\examples{

#####################################################
# Test case: the non-monotonic Sobol g-function
n <- 1000
X <- data.frame(matrix(runif(8 * n), nrow = n))
x <- johnsonshap(model = sobol.fun, X = X, N = n)
print(x)
plot(x)
library(ggplot2) ; ggplot(x)

}

\seealso{
\code{\link{johnson}}
}

\keyword{regression}
