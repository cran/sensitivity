\name{support}
\alias{support}

\title{
Support index functions: Measuring the effect of input variables over their support
}

\description{
Function to estimate the first-order and total support index functions (Fruth et al., 2016).
}

\usage{
support(model, X, Xnew = NULL, fX = NULL, gradfX = NULL, h = 1e-06, ...)
}

\arguments{
  \item{model}{a function, or a model with a predict method, defining the model to analyze.}
  \item{X}{a random sample.}
  \item{Xnew}{an optional set of points where to visualize the support indices. If missing, \code{X} is used.}
  \item{fX}{an optional vector containing the evaluations of \code{model} at \code{X}. If missing, \code{fX} is computed by evaluating \code{model} at \code{X}.}
  \item{gradfX}{an optional vector containing the evaluations of the gradient of \code{model} at \code{X}. If missing, \code{gradfX} is approximated by finite differences of \code{model} at \code{X}.}
  \item{h}{a small number for computing finite differences \code{(f(X_i + h) - f(X_i))/h}. Default is \code{1e-6}.}
  \item{\dots}{optional arguments to be passed to \code{model}.}
}

\value{
\item{main}{a matrix whose columns contain the first-order support index functions, estimated at \code{Xnew}.}
\item{total}{a matrix whose columns contain the total support index functions, estimated at \code{Xnew}.}
\item{DGSM}{a vector containing an estimation of DGSM.}
\item{X}{...}
\item{Xnew}{...}
\item{fX}{...}
\item{gradfX}{... see 'arguments' section.}
}

\details{
The first-order support index of \code{f(X)} relative to \code{X_i} is the squared conditional expectation of its partial derivative with respect to \code{X_i}.

The total support index of \code{f(X)} relative to \code{X_i} is the conditional expectation of its squared partial derivative with respect to \code{X_i}.

These two functions measure the local influence of \code{X_i}, in the global space of the other input variables. 
Up to square transformations, support indices can be viewed as regression curves of partial derivatives \code{df(X)/dX_i} with respect to \code{X_i}.
Estimation is performed by smoothing from the diagonal scatterplots \code{(X_i, df/dX_i)} with the function \code{smooth.spline{stats}} with the default options.

For the sake of comparison, support index functions may be normalized. The proposed normalization is the sum of the DGSM, equal to the sum of the overall means of total support functions. 
Normalized support index functions can be plotted with the S3 method \code{plot}, as well as the underlying diagonal scatterplots of derivatives (S3 method \code{scatterplot}).
}


\references{
J. Fruth, O. Roustant, S. Kuhnt, \emph{Support indices: Measuring the effects 
of input variables over their support}, 2016, \url{https://hal.archives-ouvertes.fr/hal-01113555}.
}

\author{
O. Roustant
}

\seealso{ 
S3 methods \code{plot} and \code{scatterplot}: \code{\link{plot.support}} 
}


\examples{

# -----------------
# ishigami function
# -----------------
n <- 5000
n.points <- 1000
d <- 3

set.seed(0)
X <- matrix(runif(d*n, min = -pi, max = pi), n, d)
Xnew <- matrix(seq(from = -pi, to = pi, length=n.points), n.points, d)

b <- support(model = ishigami.fun, X, Xnew)

# plot method (x-axis in probability scale), of the normalized support index functions
plot(b, col = c("lightskyblue4", "lightskyblue1", "black"), 
     xprob = TRUE, p = 'punif', p.arg = list(min = -pi, max = pi), ylim = c(0, 2))

# below : diagonal scatterplots of the gradient, 
# on which are based the estimation by smoothing
scatterplot(b, xprob = TRUE) 

# now with normal margins
# -----------------------
X <- matrix(rnorm(d*n), n, d)
Xnew <- matrix(rnorm(d*n.points), n.points, d)
b <- support(model = ishigami.fun, X, Xnew)

plot(b, col = c("lightskyblue4", "lightskyblue1", "black"), xprob = FALSE)
scatterplot(b, xprob = FALSE, type = "histogram", bins = 10, cex = 1, cex.lab = 1.5)
}

