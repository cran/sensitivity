\name{support}
\alias{support}
\title{
Support indices: Measuring the effect of input variables over their support
}
\description{
Function to estimate and plot the support index functions (first order and total) 
for one variable (Fruth et al., 2016) \code{xi}
}
\usage{
support(fun, d, xi = 1, h = 0.01, n = 5000, n.points = 50, q, q.arg, ...)
}

\arguments{
  \item{fun}{a function defining the model to analyze.}
  \item{d}{number of input variables of the model.}
  \item{xi}{the number of the studied input variable.}
  \item{h}{small over which to approximate the derivative.}
  \item{n}{number of Monte Carlo estimates.}
  \item{n.points}{number of equally spread points over the domain to evaluate 
    the function.}
  \item{q}{a vector of quantile functions names corresponding to wanted factors 
    distributions (see details below).}
  \item{q.arg}{a list of quantile functions parameters (see details below).}
  \item{\dots}{any other arguments for the plot of the support indices
    which are passed unchanged each time it is called.}
}

\details{
  If the arguments \code{q} and \code{q.args} are not given, the factors
  are taken uniformly distributed on \eqn{[0,1]}{[0,1]}. The
  argument \code{q} must be list of character strings, giving the names
  of the quantile functions (one for each factor), such as \code{qunif},
  \code{qnorm}\dots It can also be a single character string, meaning
  same distribution for all. The argument \code{q.arg} must be a list of
  lists, each one being additional parameters for the corresponding
  quantile function. For example, the parameters of the quantile
  function \code{qunif} could be \code{list(min=1, max=2)}, giving an
  uniform distribution on \eqn{[1,2]}{[1,2]}.  If \code{q} is a single
  character string, then \code{q.arg} must be a single list (rather than
  a list of one list).
}

\value{
  \code{support} returns a list, containing the first order support indices 
    and the total support indices
Note that the function returns plot and values at the same time.
}
\references{
J. Fruth, O. Roustant, S. Kuhnt, \emph{Support indices: Measuring the effects 
of input variables over their support}, 2016, submitted.
}
\author{
Jana Fruth
}

\examples{

col1 <- "lightskyblue1"
colT <- "lightskyblue4"
ylim <- c(0,75)
n.points <- 50  

#X1
a <- support(fun=ishigami.fun, d=3, xi=1, h=0.01, n=5000, n.points=n.points, 
        q="qunif", q.arg=list(min=-pi, max=pi), ylim=ylim, xaxt="n", ylab="$X1$")
legend("topleft", legend=c("first-order support index", "total support index"), 
  fill=c(col1,colT), bty="n")

#x2
support(fun=ishigami.fun, d=3, xi=2, h=0.01, n=5000, n.points=n.points, 
        q="qunif", q.arg=list(min=-pi, max=pi), ylim=ylim, xaxt="n", ylab="$X2$")
legend("topleft", legend=c("first-order support index", "total support index"), 
  fill=c(col1,colT), bty="n")

#x3
support(fun=ishigami.fun, d=3, xi=3, h=0.01, n=5000, n.points=n.points, 
        q="qunif", q.arg=list(min=-pi, max=pi), ylim=ylim, ylab="$X3$")
legend("topleft", legend=c("first-order support index", "total support index"), 
  fill=c(col1,colT), bty="n")
}
