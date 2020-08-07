\name{shapleyBlockEstimation}
\alias{shapleyBlockEstimation}
\alias{shapleyBlockEstimationS}
\alias{shapleyBlockEstimationX}


\title{Computation of the Shapley effects in the Gaussian linear framework
with an unknown block-diagonal covariance matrix}

\description{
\code{shapleyBlockEstimation} estimates the Shaplay effects of a Gaussian linear model
when the parameters are unknown and when the number of inputs is large,
choosing the most likely block-diagonal structure of the covariance matrix.
}

\usage{
shapleyBlockEstimationS(Beta, S, kappa=0,  M=20, tol=10^(-6))
shapleyBlockEstimationX(X, Y, delta=NULL, M=20, tol=10^(-6))
}

\arguments{
  \item{Beta}{A vector containing the (estimated) coefficients of the linear model.}
  \item{S}{Empirical covariance matrix of the inputs. Has to be positive semi-definite matrix with same size that Beta.}
  \item{X}{Matrix containing an i.i.d. sample of the inputs.}
  \item{Y}{Vector containing the corresponding i.i.d. sample of the (noisy) output.}
  \item{kappa}{The positive penalization coefficient that promotes block-diagonal matrices.}
  \item{delta}{Positive number that fixes the positive penalization coefficient
  \code{kappa} to \eqn{1/(p n^{delta})}. It is advised to choose \code{delta} to 3/2.}
  \item{M}{Maximal size of the estimate of the block-diagonal structure. The computation time grows exponentially with \code{M}.}
  \item{tol}{A relative tolerance to detect zero singular values of Sigma.}
}


\value{
  \code{shapleyBlockEstimationS} and \code{shapleyblockEstimationX} return a list containing the following compopents:
  \item{label}{a vector containing the label of the group of each input variable.}
  \item{S_B}{the block-diagonal estimated covariance matrix of the inputs.}
  \item{Shapley}{a vector containing all the estimated Shapley effects.}
}

\details{
  If \code{kappa = 0} or if \code{delta = NULL}, there is no penalization.
  
  It is advised to choose \code{M} smaller or equal than 20. For \code{M} larger or equal than 25, the computation is very long.
}

\references{
   B. Broto, F. Bachoc, L. Clouvel and J-M Martinez, 2020.\emph{Block-diagonal covariance estimation
   and application to the Shapley effects in sensitivity analysis}. arXiv preprint arXiv:1907.12780.

   B. Broto, F. Bachoc, M. Depecker, and J-M. Martinez, 2019, \emph{Sensitivity indices
for independent groups of variables}, Mathematics and Computers in Simulation, 163, 19--31.
   
   B. Iooss and C. Prieur, 2019, \emph{Shapley effects for sensitivity analysis with 
correlated inputs: comparisons with Sobol' indices, numerical estimation and 
applications}, International Journal of Uncertainty Quantification, 9, 493--514.

  A.B. Owen and C. Prieur, 2016, \emph{On Shapley value for measuring importance
of dependent inputs}, SIAM/ASA Journal of Uncertainty Quantification, 5, 986--1002.
}

\author{
Baptiste Broto, CEA LIST
}

\seealso{
\link{shapleyLinearGaussian}, \link{shapleyPermEx}, \link{shapleyPermRand}, \link{shapleySubsetMc}
}

\examples{

\donttest{
# packages for the plots of the matrices
library(ggplot2)
library(graphics)
library(gplots)


# the following function improves the plots of the matrices
sig=function(x,alpha=1)
{
  return(1/(1+exp(-x/alpha)))
}


# 1) we generate the parameters by groups in order

K=4 # number or groups

pk=rep(0,K)
for(k in 1:K)
{
  pk[k]=round(6+4*runif(1))
}
p=sum(pk)
Sigma_ord=matrix(0,nrow=p, ncol=p)
ind_min=0
L=5
for(k in 1:K){
  p_k=pk[k]
  ind=ind_min+(1:p_k)
  ind_min=ind_min+p_k
  
  A=2*matrix(runif(p_k*L),nrow=L,ncol=p_k)-1
  Sigma_ord[ind,ind]=t(A)\%*\%A + 0.2*diag(rep(1,p_k))
}


image((0:p)+0.5,(0:p)+0.5,z=sig(Sigma_ord,1),col=rev(redblue(100)), zlim=c(0,1),
      ylim=c(p+0.5,0.5), main=expression(Sigma["order"]), cex.main=3,ylab = "", 
      xlab = "",axes=FALSE)
box()


Beta_ord=3*runif(p)+1
eta_ord=shapleyLinearGaussian(Beta=Beta_ord, Sigma=Sigma_ord)
barplot(eta_ord,main=expression(eta["order"]),cex.axis = 2,cex.main=3)


# 2) We sample the input variables to get an input vector more general

samp=sample(1:p)
Sigma=Sigma_ord[samp,samp]

image((0:p)+0.5,(0:p)+0.5,z=sig(Sigma,1),col=rev(redblue(100)), zlim=c(0,1),
      ylim=c(p+0.5,0.5), main=expression(Sigma), cex.main=3,ylab = "",xlab = "",axes=FALSE)
box()


Beta=Beta_ord[samp]
eta=shapleyLinearGaussian(Beta=Beta, Sigma=Sigma)
barplot(eta,main=expression(eta),cex.axis = 2,cex.main=3)




# 3) We generate the observations with these parameters

n=5*p #sample size


C=chol(Sigma)
X0=matrix(rnorm(p*n),ncol=p)
X=X0\%*\%C

S=var(X) #empirical covariance matrix
image((0:p)+0.5,(0:p)+0.5,z=sig(S,1),col=rev(redblue(100)), zlim=c(0,1),
      ylim=c(p+0.5,0.5), main=expression(S), cex.main=3,ylab = "", xlab = "",
      axes=FALSE)
box()

beta0=rnorm(1)
Y=X\%*\%as.matrix(Beta)+beta0+0.2*rnorm(p)



# 4) We estimate the block-diagonal covariance matrix and the Shapley effects using the observations
# We assume that we know that the groups are smaller than 15

Estim=shapleyBlockEstimationX(X,Y,delta=3/4, M=15, tol=10^(-6))

eta_hat=Estim$Shapley
S_B=Estim$S_B

image((0:p)+0.5,(0:p)+0.5,z=sig(S_B,1),col=rev(redblue(100)), zlim=c(0,1),
      ylim=c(p+0.5,0.5), main=expression(S[hat(B)]), cex.main=3,ylab = "",
      xlab = "",axes=FALSE)
box()

barplot(eta_hat,main=expression(hat(eta)),cex.axis = 2,cex.main=3)


sum(abs(eta_hat-eta))
}

}
