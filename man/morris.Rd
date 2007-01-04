\name{morris}
\alias{morris}
\alias{tell.morris}
\alias{print.morris}
\alias{plot.morris}

\title{The Morris OAT Screening Method}

\description{
  \code{morris} is the implementation of the Morris OAT Screening
  method. This function generates the Morris design of experiments and
  computes the measures of sensitivity \eqn{\mu^*}{mu*} and
  \eqn{\sigma}{sigma}.
}

\usage{
morris(model = NULL, factors, levels, R, jump = NULL,
       min = 0, max = 1, scale = TRUE, optim = NULL, \dots)
}

\arguments{
  \item{model}{the model}
  \item{factors}{the number of factors, or their names}
  \item{levels}{the number of levels of the design grid}
  \item{R}{the number of repetitions of the design, i.e. the number of
    elementary effect computed per factor}
  \item{jump}{the \sQuote{grid jump} coefficient}
  \item{min}{the minimum values for the factors}
  \item{max}{the maximum values for the factors}
  \item{scale}{logical. If \code{TRUE}, the input and output data are
    scaled}
  \item{optim}{optimization of the design for better coverage of the
    space (cf Campolongo 2005), not documented yet (for informations
    feel free to ask the maintainer)}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called}
}

\details{
  \code{model} is a function or a predictor (a class with a
  \code{predict} method) computing the response \code{y} based on the
  sample given by \code{x}. If no model is specified, the indices will be
  computed when one gives the response.

  \code{factors} could either be a single number or a vector of
  character strings.

  The number of levels is not necessary the same for each space
  coordinate. It is the case when \code{levels} is a single integer.
  
  \code{min} and \code{max} are boundaries of the region of
  experimentation. They can be single values (the same for each
  factor) or vectors.

  \code{jump} is such that:
  \deqn{
    \Delta_i = \texttt{jump}_i \times \frac{\texttt{max}_i-\texttt{min}_i}{\texttt{levels}_i-1}
  }{
    Delta[i] = jump[i] * ( max[i] - min[i] ) / ( levels[i] - 1 )
  }
  If \code{jump} is given as \code{NULL} and the number of levels is
  even (for each component), then \code{jump} has the value recommended
  by Morris:
  \eqn{\texttt{jump} = \texttt{levels} / 2}{jump = levels / 2}.
  If \code{jump} is a single value, then it is taken the same for each
  coordinate.

}

\value{
  \code{morris} returns an object of class \code{"morris"}.
  An object of class \code{"morris"} is a list containing the following
  components:

  \item{x}{the design of experiments (input sample)}
  \item{y}{the response}
  \item{ee}{the matrix of the elementary effects}
  \item{mu}{the estimations of the \eqn{\mu^*}{mu*} index}
  \item{sigma}{the estimations of the \eqn{\sigma}{sigma} index}
}

\section{Computational cost}{
  The number of model evaluations is \eqn{(p+1) \times R}{(p + 1) * R}
  where \eqn{p}{p} is the number of factors.
}

\references{
Saltelli, A., Chan, K. and Scott, E. M., 2000, \emph{Sensitivity
  analysis}. Wiley.

Morris, M. D., 1991, Factorial sampling plans for preliminary
computational experiments. \emph{Technometrics}, \bold{33}, 161--174.
}

\examples{
# Test case : the non-monotonic function of Morris

sa <- morris(model = morris.fun, factors = 20, levels = 4, R = 4)
print(sa)
plot(sa)
}

\keyword{misc}
