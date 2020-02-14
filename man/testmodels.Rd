\name{testmodels}
\alias{testmodels}
\alias{sobol.fun}
\alias{ishigami.fun}
\alias{morris.fun}
\alias{atantemp.fun}
\alias{campbell1D.fun}
\alias{linkletter.fun}
\alias{heterdisc.fun}

\title{Test Models for Sensitivity Analysis}

\description{
  These functions are standard testcase for sensitivity analysis
  benchmarks. For a scalar output (see Saltelli et al. 2000, section 2.9):
  \itemize{
    \item the g-function of Sobol' with 8 inputs, X ~ U[0,1];
    \item the function  of Ishigami with 3 inputs, X ~ U[-pi,pi];
    \item the function of Morris with 20 inputs, X ~ U[0,1];
    \item the Linkletter et al. (2006) decreasing coefficients function, X ~ U[0,1];
    \item the heterdisc function with 4 inputs, X ~ U[0,20].
  }
  For functional output cases:
  \itemize{
    \item the Arctangent temporal function with 2 inputs, X ~ U[-7,7] 
    (Auder, 2011). The functional support is on [0,2pi];
    \item the Cambell1D function with 4 inputs, X ~U[-1,5] 
    (Campbell et al. 2006). The functional support is on [-90,90].
}
}

\usage{
sobol.fun(X)
ishigami.fun(X)
morris.fun(X)
atantemp.fun(X, q = 100)
campbell1D.fun(X, theta = -90:90)
linkletter.fun(X)
}

\arguments{
  \item{X}{a matrix (or \code{data.frame}) containing the input
    sample.}
    \item{q}{for the atantemp() function: the number of discretization steps of the functional output}
    \item{theta}{for the campbell1D() function: the discretization steps (angles in degrees)}
}

\value{
  A vector of function responses.
}

\references{
A. Saltelli, K. Chan and E. M. Scott eds, 2000, \emph{Sensitivity   Analysis}, Wiley.
}

\author{
Gilles Pujol and Bertrand Iooss
}

\examples{
\donttest{

# Examples for the functional toy fonctions

# atantemp function

y0 <- atantemp.fun(matrix(c(-7,0,7,-7,0,7),ncol=2))
plot(y0[1,],type="l")
apply(y0,1,lines)

n <- 100
X <- matrix(c(runif(2*n,-7,7)),ncol=2)
y <- atantemp.fun(X)
x11()
plot(y0[2,],ylim=c(-2,2),type="l")
apply(y,1,lines)

# campbell1D function

N1=100         # nombre de simulations pour courbes 1D
min=-1 ; max=5
nominal=(max+min)/2

X1 = NULL ; y1 = NULL
Xnom=matrix(nominal,nr=1,nc=4)
ynom=campbell1D.fun(Xnom,theta=-90:90)
x11()
plot(ynom,ylim=c(8,30),type="l",col="red")
for (i in 1:N1){
  X=matrix(runif(4,min=min,max=max),nr=1,nc=4)
  rbind(X1,X)
  y=campbell1D.fun(X,theta=-90:90)
  rbind(y1,y)
  lines(y)
}

}
}

\keyword{misc}
