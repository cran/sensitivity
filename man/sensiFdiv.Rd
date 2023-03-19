\name{sensiFdiv}
\alias{sensiFdiv}
\alias{tell.sensiFdiv}
\alias{print.sensiFdiv}
\alias{plot.sensiFdiv}
\alias{ggplot.sensiFdiv}

\title{Sensitivity Indices based on Csiszar f-divergence}

\description{ \code{sensiFdiv} conducts a density-based sensitivity 
analysis where the impact of an input variable is defined 
in terms of dissimilarity between the original output density function 
and the output density function when the input variable is fixed. 
The dissimilarity between density functions is measured with Csiszar f-divergences. 
Estimation is performed through kernel density estimation and 
the function \code{kde} of the package \code{ks}.}

\usage{
sensiFdiv(model = NULL, X, fdiv = "TV", nboot = 0, conf = 0.95, ...)
\method{tell}{sensiFdiv}(x, y = NULL, \dots)
\method{print}{sensiFdiv}(x, \dots)
\method{plot}{sensiFdiv}(x, ylim = c(0, 1), \dots)
\method{ggplot}{sensiFdiv}(data, mapping = aes(), ylim = c(0, 1), \dots, environment
                 = parent.frame())
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{X}{a matrix or \code{data.frame} representing the input random sample.}
  \item{fdiv}{a string or a list of strings specifying the Csiszar f-divergence 
    to be used. Available choices are "TV" (Total-Variation), 
  "KL" (Kullback-Leibler), "Hellinger" and "Chi2" (Neyman chi-squared).}
  \item{nboot}{the number of bootstrap replicates}
  \item{conf}{the confidence level for confidence intervals.}
  \item{x}{a list of class \code{"sensiFdiv"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{data}{a list of class \code{"sensiFdiv"} storing the state of the
  sensitivity study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{ylim}{y-coordinate plotting limits.}
  \item{mapping}{Default list of aesthetic mappings to use for plot. If not specified, 
    must be supplied in each layer added to the plot.}
  \item{environment}{[Deprecated] Used prior to tidy evaluation.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
}

\value{
\code{sensiFdiv} returns a list of class \code{"sensiFdiv"}, containing all
  the input arguments detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments.}
  \item{y}{a vector of model responses.}
  \item{S}{the estimations of the Csiszar f-divergence sensitivity indices. 
    If several divergences have been selected, Sis a list where each element 
    encompasses the estimations of the sensitivity indices for one of the divergence.}
}

\details{
	Some of the Csiszar f-divergences produce sensitivity indices that have 
    already been studied in the context of sensitivity analysis.
	In particular, "TV" leads to the importance measure proposed by Borgonovo (2007) 
  (up to a constant), "KL" corresponds to the mutual information (Krzykacz-Hausmann 2001) and 
	"Chi2" produces the squared-loss mutual information. See Da Veiga (2015) for details.
}
\references{
Borgonovo E. (2007), \emph{A new uncertainty importance measure}, 
Reliability Engineering and System Safety 92(6), 771--784.

Da Veiga S. (2015), \emph{Global sensitivity analysis with dependence measures},
Journal of Statistical Computation and Simulation, 85(7), 1283--1305. 

Krzykacz-Hausmann B. (2001), \emph{Epistemic sensitivity analysis based on the 
concept of entropy}, Proceedings of SAMO2001, 53--57.
}

\author{
Sebastien Da Veiga, Snecma
}

\seealso{
\code{\link[ks]{kde}, \link{sensiHSIC}}
}

\examples{
\donttest{
library(ks)

# Test case : the non-monotonic Sobol g-function
n <- 100
X <- data.frame(matrix(runif(8 * n), nrow = n))

# Density-based sensitivity analysis
# the next lines are put in comment because too long for CRAN tests
#x <- sensiFdiv(model = sobol.fun, X = X, fdiv = c("TV","KL"), nboot=30)
#print(x)
#library(ggplot2)
#ggplot(x)
}
}

