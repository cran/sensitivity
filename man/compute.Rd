\name{compute}
\alias{compute}

\title{Computation Of Sensitivity Indices}

\description{
  \code{compute} is used to compute the sensitivity indices. It is a
  generic function which means that each class of sensitivity analysis
  overloads this function. When passing a model to the sensitivity
  analysis function (such as \code{linsa}), \code{compute} is
  automatically called. \code{compute} is useful when the
  model is external (see example).
}

\usage{
compute(sa, y = NULL)
}

\arguments{
  \item{sa}{The sensitivity analysis object.}
  \item{y}{The response.}
}

\details{
  \code{sa} is an object returned by a sensitivity analysis function,
  such as \code{\link{linsa}}, \code{\link{morris}},
  \code{\link{sobol}}, \code{\link{sobol.sal02}} or \code{\link{fast}}.

  \code{y} should be a \code{numeric} vector.
}

\value{
  \code{compute} doesn't return anything, but it stores the results into
  the \code{sa} variable (it modifies its argument).
}

\seealso{
  \code{\link{sensitivity}}
}

\examples{
# Example of the FAST method

sa <- fast(model = NULL, factors = 8, n = 1000, nboot = 100)

# at this stage, only the design of experiment (sa$x) was generated
# the response is computed "manually":

y <- sobol.fun(sa$x) # at this place could be a
                     # call to an external code

# then, the sensitivity analysis:

compute(sa, y)
print(sa)
}

\keyword{misc}
