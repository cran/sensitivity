\name{PLIsuperquantile}
\alias{PLIsuperquantile}

\title{Perturbed-Law based sensitivity Indices (PLI) for superquantile}

\description{
  PLIsuperquantile computes the Perturbed-Law based Indices (PLI) for superquantile, 
  which are robustness indices related to a superquantile of a 
  model output, estimated by a Monte Carlo method. 
  See Iooss et al. (2020).
}

\usage{
PLIsuperquantile(order,x,y,deltasvector,InputDistributions,type="MOY",
  samedelta=TRUE, percentage=TRUE,nboot=0,conf=0.95,bootsample=TRUE,bias=TRUE)
}

\arguments{
  \item{order}{the order of the superquantile to estimate.}
  \item{x}{the matrix of simulation points coordinates, 
    one column per variable.}
    \item{y}{the vector of model outputs.}
	\item{deltasvector}{a vector containing the values of delta for which
    the indices will be computed.}
	\item{InputDistributions}{a list of list. Each list contains, as a list,
    the name of the distribution to be used and the parameters.
	  Implemented cases so far:
	  \itemize{
      \item For a mean perturbation: Gaussian, Uniform, Triangle,
        Left Trucated Gaussian, Left Truncated Gumbel. Using Gumbel
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
  \item{percentage}{a boolean that defines the formula used for the PLI.
    \itemize{
      \item If it is set at FALSE, the classical formula used in the bibliographical 
      references is used.
      \item If not (set as TRUE), the PLI is given in percentage of variation 
      of the superquantile (even if it is negative).
        }}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level for bootstrap confidence intervals.}
  \item{bootsample}{If TRUE, the uncertainty about the original quantile estimation
      is taken into account in the PLI confidence intervals (see Iooss et al., 2020). 
      If FALSE, standard confidence intervals are computed for the PLI.
      It mainly changes the CI at small delta values.}
  \item{bias}{defines the version of PLI-superquantile:
    \itemize{
      \item If it is set at "TRUE", it gives the mean of outputs above the perturbed 
      quantile (alternative formula)
      \item If it is set at "FALSE", it gives the mean of perturbed outputs above the 
      perturbed quantile (original formula)
        }}
}

\value{
  \code{PLIsuperquantile} returns a list of matrices (each column corresponds to an input, 
  each line corresponds to a twist of amplitude delta) 
  containing the following components:
  \item{PLI}{the PLI.}
  \item{PLICIinf}{the bootstrap lower confidence interval values of the PLI.}
  \item{PLICIsup}{the bootstrap upper confidence interval values of the PLI.}
  \item{superquantile}{the perturbed superquantile.}
  \item{superquantileCIinf}{the bootstrap lower confidence interval values of the perturbed superquantile.}
  \item{superquantileCIsup}{the bootstrap upper confidence interval values of the perturbed superquantile.}
}

\references{
B. Iooss, V. Verges and V. Larget, 2022, \emph{BEPU robustness analysis via perturbed law-based 
sensitivity indices}, Proceedings of the Institution of Mechanical Engineers, 
Part O: Journal of Risk and Reliability, 236:855-865.

P. Lemaitre, E. Sergienko, A. Arnaud, N. Bousquet, F. Gamboa and B. Iooss, 2015, 
\emph{Density modification based reliability sensitivity analysis}, Journal of Statistical 
Computation and Simulation, 85:1200-1223. 

}

\author{
Bertrand Iooss
}

\seealso{
  \code{\link{PLI}, \link{PLIquantile}, \link{PLIsuperquantile_multivar}}
}

\examples{
\donttest{

# Model: 3D function 

  distribution = list()
	for (i in 1:3) distribution[[i]]=list("norm",c(0,1))
  
# Monte Carlo sampling 

  N = 10000
	X = matrix(0,ncol=3,nrow=N)
	for(i in 1:3) X[,i] = rnorm(N,0,1)
     
	Y = 2 * X[,1] + X[,2] + X[,3]/2
	alpha <- 0.95
	
	q95 = quantile(Y,alpha)
  sq95a <- mean(Y*(Y>q95)/(1-alpha)) ; sq95b <- mean(Y[Y>q95])
	
	nboot=20 # change to nboot=200 for consistency
	
# sensitivity indices with perturbation of the mean 
  
	v_delta = seq(-1,1,1/10) 
	toto = PLIsuperquantile(alpha,X,Y,deltasvector=v_delta,
	  InputDistributions=distribution,type="MOY",samedelta=TRUE,
	  percentage=FALSE,nboot=nboot,bias=TRUE)

# Plotting the PLI
  par(mar=c(4,5,1,1))
	plot(v_delta,toto$PLI[,2],ylim=c(-0.5,0.5),xlab=expression(delta),
		ylab=expression(hat(PLI[i*delta])),pch=19,cex=1.5)
	points(v_delta,toto$PLI[,1],col="darkgreen",pch=15,cex=1.5)
	points(v_delta,toto$PLI[,3],col="red",pch=17,cex=1.5)
	lines(v_delta,toto$PLICIinf[,2],col="black")
	lines(v_delta,toto$PLICIsup[,2],col="black")
	lines(v_delta,toto$PLICIinf[,1],col="darkgreen")
	lines(v_delta,toto$PLICIsup[,1],col="darkgreen")
	lines(v_delta,toto$PLICIinf[,3],col="red")
	lines(v_delta,toto$PLICIsup[,3],col="red")
	abline(h=0,lty=2)
	legend(-1,0.5,legend=c("X1","X2","X3"),
		col=c("darkgreen","black","red"),pch=c(15,19,17),cex=1.5)
  
# Plotting the perturbed superquantiles
  par(mar=c(4,5,1,1))
	plot(v_delta,toto$superquantile[,2],ylim=c(3,7),xlab=expression(delta),
		ylab=expression(hat(q[i*delta])),pch=19,cex=1.5)
	points(v_delta,toto$superquantile[,1],col="darkgreen",pch=15,cex=1.5)
	points(v_delta,toto$superquantile[,3],col="red",pch=17,cex=1.5)
	lines(v_delta,toto$superquantileCIinf[,2],col="black")
	lines(v_delta,toto$superquantileCIsup[,2],col="black")
	lines(v_delta,toto$superquantileCIinf[,1],col="darkgreen")
	lines(v_delta,toto$superquantileCIsup[,1],col="darkgreen")
	lines(v_delta,toto$superquantileCIinf[,3],col="red")
	lines(v_delta,toto$superquantileCIsup[,3],col="red")
	abline(h=q95,lty=2)
	legend(-1,7,legend=c("X1","X2","X3"),
		col=c("darkgreen","black","red"),pch=c(15,19,17),cex=1.5)
		
# Plotting the unbiased PLI in percentage with refined confidence intervals
	toto = PLIsuperquantile(alpha,X,Y,deltasvector=v_delta,
	  InputDistributions=distribution,type="MOY",samedelta=TRUE,percentage=TRUE,
	  nboot=nboot,bootsample=FALSE,bias=FALSE)
	  
  par(mar=c(4,5,1,1))
	plot(v_delta,toto$PLI[,2],ylim=c(-0.4,0.5),xlab=expression(delta),
		ylab=expression(hat(PLI[i*delta])),pch=19,cex=1.5)
	points(v_delta,toto$PLI[,1],col="darkgreen",pch=15,cex=1.5)
	points(v_delta,toto$PLI[,3],col="red",pch=17,cex=1.5)
	lines(v_delta,toto$PLICIinf[,2],col="black")
	lines(v_delta,toto$PLICIsup[,2],col="black")
	lines(v_delta,toto$PLICIinf[,1],col="darkgreen")
	lines(v_delta,toto$PLICIsup[,1],col="darkgreen") 
	lines(v_delta,toto$PLICIinf[,3],col="red")
	lines(v_delta,toto$PLICIsup[,3],col="red")
	abline(h=0,lty=2)
	legend(-1,0.5,legend=c("X1","X2","X3"),
		col=c("darkgreen","black","red"),pch=c(15,19,17),cex=1.5)

##################################################
# another visualization by using the plotCI() fct 
# (from plotrix package) for the CI plotting (from Vanessa Verges)

	library(plotrix)
	parameters = list(colors=c("darkgreen","black","red"),symbols=c(15,19,17),
	  overlay=c(FALSE,TRUE,TRUE))
  par(mar=c(4,5,1,1),xpd=TRUE)
  for (i in 1:3){
  plotCI(v_delta,toto$PLI[,i],ui=toto$PLICIsup[,i],li=toto$PLICIinf[,i],
         cex=1.5,col=parameters$colors[i],pch=parameters$symbols[i],
         add=parameters$overlay[i], xlab="", ylab="")
  }
  title(xlab=expression(delta),ylab=expression(hat(PLI[i*delta])),
      main=bquote("PLI-superquantile (N ="~.(N) ~ ","~alpha~"="~.(alpha)~
      ") of Y="~2*X[1] + X[2] + X[3]/2))
  abline(h=0,lty=2)
  legend("topleft",legend=c("X1","X2","X3"),
          col=parameters$colors,pch=parameters$symbols,cex=1.5)

}
}
