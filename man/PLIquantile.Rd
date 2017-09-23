\name{PLIquantile}
\alias{PLIquantile}

\title{Perturbed-Law based sensitivity Indices (PLI) for quantile}

\description{
  PLIquantile computes the Perturbed-Law based Indices (PLI) for quantile, 
  which are sensitivity indices related to a quantile of a 
  model output, estimated by a Monte Carlo method, See Sueur et al. (2016, 2017).
}

\usage{
PLIquantile(order,x,y,quantilehat,deltasvector,
       InputDistributions,type="MOY",samedelta=TRUE)
}

\arguments{
  \item{order}{the order of the quantile to estimate.}
  \item{x}{the matrix of simulation points coordinates, 
    one column per variable.}
    \item{y}{the vector of model outputs.}
  \item{quantilehat}{the estimation of quantile q.}
	\item{deltasvector}{a vector containing the values of delta for which
    the indices will be computed.}
	\item{InputDistributions}{a list of list. Each list contains, as a list,
    the name of the distribution to be used and the parameters.
	  Implemented cases so far:
	  \itemize{
      \item For a mean perturbation: Gaussian, Uniform, Triangle,
        Left Trucated Gaussian, Left Truncated Gumbel. Using Gumber
        requires the package \code{evd}.
	    \item For a variance perturbation: Gaussian, Uniform.
      }}
	\item{type}{a character string in which the user will specify the type of 
    perturbation wanted. 
	  The sense of "deltasvector" varies according to the type of perturbation:
    \itemize{
      \item type can take the value "MOY",in which case deltasvector is a 
        vector of perturbated means.
	    \item type can take the value "VAR",in which case deltasvector is a 
        vector of perturbated variances, therefore needs to be positive integers.
        }}
	\item{samedelta}{a boolean used with the value "MOY" for type. 
    \itemize{
      \item If it is set at TRUE, the mean perturbation will be the same for all  
        the variables. 
      \item If not, the mean perturbation will be new_mean = mean+sigma*delta 
        where mean, sigma are parameters defined in InputDistributions and 
        delta is a value of deltasvector.
        }}
}

\value{
  \code{PLIquantile} returns a matrix where the PLI are stored. Each column corresponds to 
      an input, each line corresponds to a twist of amplitude delta.
}

\references{
P. Lemaitre, E. Sergienko, A. Arnaud, N. Bousquet, F. Gamboa and B. Iooss, 2015, 
\emph{Density modification based reliability sensitivity analysis}, Journal of Statistical 
Computation and Simulation, 85:1200-1223. 

R. Sueur, N. Bousquet, B. Iooss and J. Bect, 2016,
\emph{Perturbed-Law based sensitivity Indices for sensitivity analysis in structural reliability},
Proceedings of the SAMO 2016 Conference, Reunion Island, France, December 2016.

R. Sueur, B. Iooss and T. Delage, 2017,
\emph{Sensitivity analysis using perturbed-law based indices for quantiles and application to an industrial case}, 
10th International Conference on Mathematical Methods in Reliability (MMR 2017), Grenoble, France, July 2017.
}

\author{
Paul Lemaitre, Bertrand Iooss and Thibault Delage
}

\seealso{
  \code{\link{PLI}}
}

\examples{
\dontrun{

# Model: 3D function 

  distribution = list()
	for (i in 1:3) distribution[[i]]=list("norm",c(0,1))
  
# Monte Carlo sampling to obtain failure points

  N = 10000
	X = matrix(0,ncol=3,nrow=N)
	for(i in 1:3) X[,i] = rnorm(N,0,1)
     
	Y = 2 * X[,1] + X[,2] + X[,3]/2
	q95 = quantile(Y,0.95)
	
# sensitivity indices with perturbation of the mean 
  
	v_delta = seq(-1,1,1/10) 
	toto = PLIquantile(0.95,X,Y,q95,deltasvector=v_delta,
	  InputDistributions=distribution,type="MOY",samedelta=TRUE)

	par(mar=c(4,5,1,1))
	plot(v_delta,toto[,2],ylim=c(-4.5,4.5),xlab=expression(delta),
		ylab=expression(hat(S[i*delta])),pch=19,cex=1.5)
	points(v_delta,toto[,1],col="darkgreen",pch=15,cex=1.5)
	points(v_delta,toto[,3],col="red",pch=17,cex=1.5)
	abline(h=0,lty=2)
	legend(0.8,4.4,legend=c("X1","X2","X3"),
		col=c("darkgreen","black","red"),pch=c(15,19,17),cex=1.5)
  
}
}
