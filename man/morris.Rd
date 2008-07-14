\name{morris}
\alias{morris}
\alias{tell.morris}
\alias{print.morris}
\alias{plot.morris}
\alias{plot3d.morris}

\title{Morris's Elementary Effects Screening Method}

\description{
  
  \code{morris} implements the Morris's elementary effects screening
  method (Morris 1992). This method, based on design of experiments,
  allows to identify the few important factors at a cost of \eqn{r
  \times (p+1)}{r * (p + 1)} simulations (where \eqn{p}{p} is the number
  of factors). This implementation includes some improvements of the
  original method: space-filling optimization of the design (Campolongo
  et al. 2007) and simplex-based design (Pujol 2008).
}

\usage{
morris(model = NULL, factors, r, design, binf = 0, bsup = 1,
       scale = TRUE, \dots)
\method{tell}{morris}(x, y = NULL, \dots)
\method{print}{morris}(x, \dots)
\method{plot}{morris}(x, identify = FALSE, \dots)
plot3d.morris(x, alpha = c(0.2, 0), sphere.size = 1)
}

\arguments{
  \item{model}{a function, or a model with a \code{predict} method,
    defining the model to analyze.}
  \item{factors}{an integer giving the number of factors, or a vector of
    character strings giving their names.}
  \item{r}{either an integer giving the number of repetitions of the design,
    i.e. the number of elementary effect computed per factor, or a
    vector of two integers \code{c(r1, r2)} for the space-filling
    improvement (Campolongo et al.). In this case, \code{r1} is the
    wanted design size, and \code{r2} (\eqn{> \code{r1}}{> \code{r1}})
    is the size of the (bigger) population in which is extracted the
    design (this can throw a warning, see below).
  }
  \item{design}{a list specifying the design type and its
    parameters:
    \itemize{
      \item \code{type = "oat"} for Morris's OAT design (Morris 1992),
      with the parameters:
      \itemize{
	\item \code{levels} : either an integer specifying the number of
	levels of the design, or a vector of integers for different
	values for each factor.
	\item \code{grid.jump} : either an integer specifying the number of
	levels that are increased/decreased for computing the
	elementary effects, or a vector of integers for different values
	for each factor. If not given, it is set to \code{grid.jump =
	  1}. Notice that this default value of one does not follow
	Morris's recommendation of 
	\eqn{\texttt{levels} / 2}{levels / 2}.
      }
      \item \code{type = "simplex"} for simplex-based design (Pujol
      2008), with the parameter:
      \itemize{
	\item \code{scale.factor} : a numeric value, the homothety factor of
	the (isometric) simplexes. Edges equal one with a scale factor
	of one.
      }
    }
  }
  \item{binf}{either an integer, specifying the minimum value for the
    factors, or a vector for different values for each factor.}
  \item{bsup}{either an integer, specifying the maximum value for the
    factors, or a vector for different values for each factor.}
  \item{scale}{logical. If \code{TRUE}, the input and output data are
    scaled before computing the elementary effects.}
  \item{x}{a list of class \code{"morris"} storing the state of the
    screening study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{identify}{logical. If \code{TRUE}, the user selects with the
    mouse the factors to label on the \eqn{(\mu^*,\sigma)}{(mu*, sigma)}
  graph (see \code{identify}).}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
  \item{alpha}{a vector of three values between 0.0 (fully transparent) and 1.0
    (opaque) (see \code{rgl.material}). The first value is for the
    cone, the second for the planes.}
  \item{sphere.size}{a numeric value, the scale factor for displaying the
    spheres.}
}

\details{

  \code{plot2d} draws the \eqn{(\mu^*,\sigma)}{(mu*, sigma)} graph.
  
  \code{plot3d.morris} draws the \eqn{(\mu, \mu^*,\sigma)}{(mu, mu*,
    sigma)} graph (requires the \pkg{rgl} package). On this graph, the
    points are in a domain bounded by a cone and two planes (application
    of the Cauchy-Schwarz inequality).
  }

\value{
  \code{morris} returns a list of class \code{"morris"}, containing all
  the input argument detailed before, plus the following components:

  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments.}
  \item{y}{a vector of model responses.}
  \item{ee}{a \eqn{r \times p}{r * p} matrix of elementary effects for
    all the factors.}

  Notice that the statitics of interest (\eqn{\mu}{mu}, \eqn{\mu^*}{mu*}
  and \eqn{\sigma}{sigma}) are not stored. They can be printed by the
  \code{print} method, but to extract numerical values, one has to
  compute them with the following instructions:
  \preformatted{
mu <- apply(x$ee, 2, mean)
mu.star <- apply(x$ee, 2, function(x) mean(abs(x)))
sigma <- apply(x$ee, 2, sd)}
}

\section{Warning messages}{
  \describe{
    \item{"keeping r' repetitions out of r"}{when generating the design of
  experiments, identical repetitions are removed, leading to a lower
  number than requested.}
  }

}

\references{
M. D. Morris, 1991, \emph{Factorial sampling plans for preliminary
computational experiments}, Technometrics, 33, 161--174.

F. Campolongo, J. Cariboni and A. Saltelli, 2007, \emph{An effective
screening design for sensitivity}, Environmental Modelling \&
  Software, 22, 1509--1518.

G. Pujol (2008), \emph{Simplex-based screening designs for estimating
  metamodels}, submited to Reliability Engineering and System Safety.
}

\seealso{
  \code{\link{decoupling}}
}

\examples{
# Test case : the non-monotonic function of Morris
x <- morris(model = morris.fun, factors = 20, r = 4,
            design = list(type = "oat", levels = 5, grid.jump = 3))
print(x)
plot(x)
\dontrun{morris.plot3d(x)  # (requires the package 'rgl')}
}

\keyword{design}
