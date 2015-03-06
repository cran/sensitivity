\name{PLI}
\alias{PLI}

\title{Perturbation Law-based Indices (PLI)}

\description{
  PLI computes the Perturbation Law-based Indices (PLI), also known as the
  Density Modification Based Reliability Sensitivity Indices (DMBRSI), 
  which are sensitivity indices related to a probability of exceedence of a 
  model output (i.e. a failure probability), estimated by a Monte Carlo method. 
  See Lemaitre et al. (2015).
}

\usage{
PLI(failurepoints,failureprobabilityhat,samplesize,deltasvector,
       InputDistributions,type="MOY",samedelta=TRUE)
}

\arguments{
  \item{failurepoints}{a matrix of failure points coordinates, 
    one column per variable.}
  \item{failureprobabilityhat}{the estimation of failure probability P 
    through rough Monte Carlo method.}
	\item{samplesize}{the size of the sample used to estimate P. 
    One must have Pchap=dim(failurepoints)[1]/samplesize}
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
  \code{PLI} returns a list of size 2, including:
  \itemize{
    \item A matrix where the PLI are stored. Each column corresponds to 
      an input, each line corresponds to a twist of amplitude delta.
    \item A matrix where their standard deviation are stored.
  }
}

\references{P. Lemaitre, E. Sergienko, A. Arnaud, N. Bousquet, F. Gamboa and B. Iooss.
\emph{Density modification based reliability sensitivity analysis}, Journal of Statistical 
Computation and Simulation, 85:1200-1223, 2015. 
\url{http://hal.archives-ouvertes.fr/hal-00737978}
}

\author{
Paul Lemaitre
}

\examples{
# Model: Ishigami function with a treshold at -7
# Failure points are those < -7

  distributionIshigami = list()
	for (i in 1:3){
		distributionIshigami[[i]]=list("unif",c(-pi,pi))
		distributionIshigami[[i]]$r=("runif")
	}
  
# Monte Carlo sampling to obtain failure points

  N = 10^5
	X = matrix(0,ncol=3,nrow=N)
	for( i in 1:3){
    X[,i] = runif(N,-pi,pi)
    }
     
	T = ishigami.fun(X)
	s = sum(as.numeric(T < -7)) # Number of failure
	pdefchap = s/N      # Failure probability
	ptsdef = X[T < -7,]	# Failure points
	
# sensitivity indices with perturbation of the mean 
  
	v_delta = seq(-3,3,1/20) 
	Toto = PLI(failurepoints=ptsdef,failureprobabilityhat=pdefchap,samplesize=N,
		deltasvector=v_delta,InputDistributions=distributionIshigami,type="MOY",
		samedelta=TRUE)
	BIshm = Toto[[1]]
	SIshm = Toto[[2]]

	par(mar=c(4,5,1,1))
	plot(v_delta,BIshm[,2],ylim=c(-4,4),xlab=expression(delta),
		ylab=expression(hat(S[i*delta])),pch=19,cex=1.5)
	points(v_delta,BIshm[,1],col="darkgreen",pch=15,cex=1.5)
	points(v_delta,BIshm[,3],col="red",pch=17,cex=1.5)
	lines(v_delta,BIshm[,2]+1.96*SIshm[,2],col="black"); 
	lines(v_delta,BIshm[,2]-1.96*SIshm[,2],col="black")
	lines(v_delta,BIshm[,1]+1.96*SIshm[,1],col="darkgreen"); 
	lines(v_delta,BIshm[,1]-1.96*SIshm[,1],col="darkgreen")
	lines(v_delta,BIshm[,3]+1.96*SIshm[,3],col="red"); 
	lines(v_delta,BIshm[,3]-1.96*SIshm[,3],col="red");
	abline(h=0,lty=2)
	legend(0,3,legend=c("X1","X2","X3"),
		col=c("darkgreen","black","red"),pch=c(15,19,17),cex=1.5)
  
# sensitivity indices with perturbation of the variance 

	v_delta = seq(1,5,1/4) # user parameter. (the true variance is 3.29)	
	Toto = PLI(failurepoints=ptsdef,failureprobabilityhat=pdefchap,samplesize=N,
		deltasvector=v_delta,InputDistributions=distributionIshigami,type="VAR",
		samedelta=TRUE)
	BIshv=Toto[[1]]
	SIshv=Toto[[2]]

	par(mfrow=c(2,1),mar=c(1,5,1,1)+0.1)
	plot(v_delta,BIshv[,2],ylim=c(-.5,.5),xlab=expression(V_f),
		ylab=expression(hat(S[i*delta])),pch=19,cex=1.5)
	points(v_delta,BIshv[,1],col="darkgreen",pch=15,cex=1.5)
	points(v_delta,BIshv[,3],col="red",pch=17,cex=1.5)
	lines(v_delta,BIshv[,2]+1.96*SIshv[,2],col="black"); 
	lines(v_delta,BIshv[,2]-1.96*SIshv[,2],col="black")
	lines(v_delta,BIshv[,1]+1.96*SIshv[,1],col="darkgreen"); 
	lines(v_delta,BIshv[,1]-1.96*SIshv[,1],col="darkgreen")
	lines(v_delta,BIshv[,3]+1.96*SIshv[,3],col="red"); 
	lines(v_delta,BIshv[,3]-1.96*SIshv[,3],col="red");

	par(mar=c(4,5.1,1.1,1.1))
	plot(v_delta,BIshv[,2],ylim=c(-30,.7),xlab=expression(V[f]),
		ylab=expression(hat(S[i*delta])),pch=19,cex=1.5)
	points(v_delta,BIshv[,1],col="darkgreen",pch=15,cex=1.5)
	points(v_delta,BIshv[,3],col="red",pch=17,cex=1.5)
	lines(v_delta,BIshv[,2]+1.96*SIshv[,2],col="black"); 
	lines(v_delta,BIshv[,2]-1.96*SIshv[,2],col="black")
	lines(v_delta,BIshv[,1]+1.96*SIshv[,1],col="darkgreen"); 
	lines(v_delta,BIshv[,1]-1.96*SIshv[,1],col="darkgreen")
	lines(v_delta,BIshv[,3]+1.96*SIshv[,3],col="red"); 
	lines(v_delta,BIshv[,3]-1.96*SIshv[,3],col="red");
	legend(2.5,-10,legend=c("X1","X2","X3"),col=c("darkgreen","black","red"),
		pch=c(15,19,17),cex=1.5)
  
}
