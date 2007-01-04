\name{tell}
\alias{tell}

\title{Computation Of Sensitivity Indices When The Model Is External}

\description{
  \code{tell} is used to tell a sensitivity analysis object the results
  of the simulations. It is used when the model is not given when
  parametring the sensitivity analysis (whith the model=NULL
  argument). For example, it is the case when the model is external to
  \R.
}

\usage{
tell(sa, y = NULL)
}

\arguments{
  \item{sa}{the sensitivity analysis object}
  \item{y}{the response}
}

\details{
  \code{sa} is an object returned by a sensitivity analysis function,
  such as \code{\link{srcpcc}}, \code{\link{morris}},
  \code{\link{sobol}},\dots

  \code{y} should be a \code{numeric} vector.
}

\value{
  \code{tell} doesn't return anything. It does the sensitivity
  analysis, and stores all the results in the variable \code{sa}.
}

\examples{
# Example of the FAST method
# (one should note the call with model = NULL)

sa <- fast(model = NULL, factors = 8, n = 1000)

# at this stage, only the design of experiment (sa$x) was generated
# the response is computed "manually":

y <- sobol.fun(sa$x) # at this place could be a
                     # call to an external code

# then, the sensitivity analysis:

tell(sa, y)
print(sa)

# Remark : because the model is a simple R function,
# this example is equivalent to :
\dontrun{sa <- fast(model = sobol.fun, factors = 8, n = 1000)}
}

\keyword{misc}
