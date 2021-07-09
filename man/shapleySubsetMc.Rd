\name{shapleySubsetMc}
\alias{shapleySubsetMc}
\alias{plot.shapleySubsetMc}

\title{Estimation of Shapley effects from data using nearest neighbors method}
  
\description{
  \code{shapleySubsetMc} implements the estimation of
  the Shapley effects from data using some nearest neighbors method
  to generate according to the conditional distributions of the inputs.
  It can be used with categorical inputs.
}

\usage{
shapleySubsetMc(X,Y, Ntot=NULL, Ni=3, cat=NULL, weight=NULL, discrete=NULL, noise=FALSE)
\method{plot}{shapleySubsetMc}(x, ylim = c(0, 1), \dots)
}

\arguments{
  \item{X}{a matrix or a dataframe of the input sample}
  \item{Y}{a vector of the output sample}
  \item{Ntot}{an integer of the approximate cost wanted}
  \item{Ni}{the number of nearest neighbours taken for each point}
  \item{cat}{a vector giving the indices of the input categorical variables}
  \item{weight}{a vector with the same length of \code{cat} giving the weight of each
     categorical variable in the product distance}
  \item{discrete}{a vector giving the indices of the input variable that are
     real, and not categorical, but that can take several times the same values}
  \item{noise}{logical. If FALSE (the default), the variable Y is a function of X}
  \item{x}{a list of class \code{"shapleySubsetMc"} storing the state of the
  sensitivity study (Shapley effects, cost, names of inputs)}
  \item{ylim}{y-coordinate plotting limits}
  \item{\dots}{any other arguments for plotting}
}

\value{
  \code{shapleySubsetMc} returns a list of class \code{"shapleySubsetMc"}, 
  containing:
  \item{shapley}{the Shapley effects estimates.}
  \item{cost}{the real total cost of these estimates: the total number of points for which
   the nearest neighbours were computed.}
  \item{names}{the labels of the input variables.}
    }

\details{
  If \code{weight = NULL}, all the categorical variables will have the same weight 1.
  
  If \code{Ntot = NULL}, the nearest neighbours will be compute for all the \eqn{n (2^p-2)} points,
  where n is the length of the sample. The estimation can be very long with this parameter. 
}

\references{
  B. Broto, F. Bachoc, M. Depecker, 2020, \emph{Variance reduction for estimation of 
  Shapley effects and adaptation to unknown input distribution}, 
  SIAM/ASA Journal of Uncertainty Quantification, 8:693-716.
}

\author{
Baptiste Broto
}

\seealso{
\link{shapleyPermEx}, \link{shapleyPermRand}, \link{shapleyLinearGaussian}, \link{sobolrank}, \link{shapleysobol_knn}
}

\examples{
\donttest{

# First example: the linear Gaussian framework

# we generate a covariance matrice Sigma
p <- 4 #dimension
A <- matrix(rnorm(p^2),nrow=p,ncol=p)
Sigma <- t(A)\%*\%A # it means t(A)%*%A
C <- chol(Sigma)
n <- 500 #sample size (put n=2000 for more consistency)

Z=matrix(rnorm(p*n),nrow=n,ncol=p)
X=Z\%*\%C # X is a gaussian vector with zero mean and covariance Sigma
Y=rowSums(X) 
Shap=shapleySubsetMc(X=X,Y=Y,Ntot=5000)
plot(Shap)


#Second example: The Sobol model with heterogeneous inputs

p=8 #dimension
A=matrix(rnorm(p^2),nrow=p,ncol=p)
Sigma=t(A)\%*\%A
C=chol(Sigma)
n=500 #sample size (put n=5000 for more consistency)

Z=matrix(rnorm(p*n),nrow=n,ncol=p)
X=Z%*%C+1 # X is a gaussian vector with mean (1,1,..,1) and covariance Sigma

#we create discrete and categorical variables
X[,1]=round(X[,1]/2) 
X[,2]=X[,2]>2
X[,4]=-2*round(X[,4])+4
X[(X[,6]>0 &X[,6]<1),6]=1

cat=c(1,2)  # we choose to take X1 and X2 as categorical variables (with the discrete distance)
discrete=c(4,6) # we indicate that X4 and X6 can take several times the same value

Y=sobol.fun(X)
Ntot <- 2000 # put Ntot=20000 for more consistency
Shap=shapleySubsetMc(X=X,Y=Y, cat=cat, discrete=discrete, Ntot=Ntot, Ni=10)

plot(Shap)
}

}