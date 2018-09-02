%Documentation for morrisMultOut.
\name{morrisMultOut}
\alias{morrisMultOut}
\alias{tell.morrisMultOut}

\title{Morris's Elementary Effects Screening Method for Multidimensional Outputs}

\description{
  
  \code{morrisMultOut} extend the Morris's elementary effects screening
  method (Morris 1991) to model with multidimensional outputs. 
}

\usage{
morrisMultOut(model = NULL, factors, r, design, binf = 0, bsup = 1,
       scale = TRUE, \dots)
\method{tell}{morrisMultOut}(x, y = NULL, \dots)
}

\arguments{
  \item{model}{NULL or a function returning a outputs a matrix having as columns the model outputs.}
  \item{factors}{an integer giving the number of factors, or a vector of
    character strings giving their names.}
  \item{r}{either an integer giving the number of repetitions of the design,
    i.e. the number of elementary effect computed per factor, or a
    vector of two integers \code{c(r1, r2)} for the space-filling
    improvement (Campolongo et al. 2007). In this case, \code{r1} is the
    wanted design size, and \code{r2} (\eqn{> \code{r1}}{> \code{r1}})
    is the size of the (bigger) population in which is extracted the
    design (this can throw a warning, see below).
  }
  \item{design}{a list specifying the design type and its
    parameters:
    \itemize{
      \item \code{type = "oat"} for Morris's OAT design (Morris 1991),
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
      2009), with the parameter:
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
  \item{scale}{logical. If \code{TRUE}, the input design of experiments is
    scaled after building the design and before computing the elementary 
    effects so that all factors vary within the range [0,1]. For each factor, 
    the scaling is done relatively to its corresponding bsup and binf.}
  \item{x}{a list of class \code{"morris"} storing the state of the
    screening study (parameters, data, estimates).}
  \item{y}{a vector of model responses.}
  \item{\dots}{for \code{morrisMultOut}: any other arguments for \code{model} which 
    are passed unchanged each time it is called. For \code{plot.morris}: 
    arguments to be passed to \code{plot.default}.}
}

\details{
  All the methods available for object of class \code{"morris"} are available also for objects of class \code{"morrisMultOut"}.
  See the documentation relative to the function \code{"morris"} for more details.
  }
\value{
  \code{morrisMultOut} returns a list of class \code{"c(morrisMultOut, morris)"}, containing all
  the input argument detailed before, plus the following components:
  \item{call}{the matched call.}
  \item{X}{a \code{data.frame} containing the design of experiments.}
  \item{y}{a matrix having as columns the model responses.}
  \item{ee}{a vector of aggregated elementary effects.}
}

\references{
  Monari F. and P. Strachan, 2017. \emph{Characterization of an airflow network model by sensitivity analysis: parameter screening, fixing, prioritizing and mapping}.
  Journal of Building Performance Simulation, 2017, 10, 17-36.
}

\author{
  Filippo Monari
}

\seealso{
  \code{\link{morris}}
}

\examples{

  mdl <- function (X) t(atantemp.fun(X))

  x = morrisMultOut(model = mdl, factors = 4, r = 50, 
  design = list(type = "oat", levels = 5, grid.jump = 3), binf = -1, bsup = 5, scale = FALSE)
  print(x)
  plot(x)

  x = morrisMultOut(model = NULL, factors = 4, r = 50, 
  design = list(type = "oat", levels = 5, grid.jump = 3), binf = -1, bsup = 5, scale = FALSE)
  Y = mdl(x[['X']])
  tell(x, Y)	
  print(x)
  plot(x)	
}
