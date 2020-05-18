\name{shapleyLinearGaussian}
\alias{shapleyLinearGaussian}


\title{Computation of the Shapley effects in the linear Gaussian framework}
  
\description{
  \code{shapleyLinearGaussian} implements the computation of
  the Shapley effects in the linear Gaussian framework, using the linear model
  (without the value at zero) and the covariance matrix of the inputs.
}

\usage{
shapleyLinearGaussian(Beta,Sigma,tol=10^(-6))
}

\arguments{
  \item{Beta}{a vector containing the coefficients of the linear model (without the value at zero).}
  \item{Sigma}{covariance matrix of the inputs. Has to be positive semi-definite matrix with same size that Beta.}
  \item{tol}{a relative tolerance to detect zero singular values of Sigma.}
}

\value{
  \code{shapleyLinearGaussian} returns a numeric vector containing all the Shapley effects.
}

\references{
   B. Broto, F. Bachoc, M. Depecker, and J-M. Martinez, 2019, \emph{Sensitivity indices
for independent groups of variables}, Mathematics and Computers in Simulation, 163, 19--31.
   
   B. Iooss and C. Prieur, 2019, \emph{Shapley effects for sensitivity analysis with 
correlated inputs: comparisons with Sobol' indices, numerical estimation and 
applications}, International Journal for Uncertainty Quantification, 9, 493--514.

  A.B. Owen and C. Prieur, 2016, \emph{On Shapley value for measuring importance
of dependent inputs}, SIAM/ASA Journal of Uncertainty Quantification, 5, 986--1002.
}

\author{
Baptiste Broto
}

\seealso{
\link{shapleyPermEx}, \link{shapleyPermRand}, \link{shapleySubsetMc}
}


\examples{

library(MASS)
library(igraph)

# First example:

p=5 #dimension
A=matrix(rnorm(p^2),nrow=p,ncol=p)
Sigma=t(A)\%*\%A
Beta=runif(p)
Shapley=shapleyLinearGaussian(Beta,Sigma)
plot(Shapley)


# Second Example, block-diagonal:

K=5 #number of groups
m=5 # number of variables in each group
p=K*m
Sigma=matrix(0,ncol=p,nrow=p)

for(k in 1:K)
{
  A=matrix(rnorm(m^2),nrow=m,ncol=m)
  Sigma[(m*(k-1)+1):(m*k),(m*(k-1)+1):(m*k)]=t(A)\%*\%A
}
# we mix the variables:
samp=sample(1:p,p)
Sigma=Sigma[samp,samp]

Beta=runif(p)
Shapley=shapleyLinearGaussian(Beta,Sigma)
plot(Shapley)

}