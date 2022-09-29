\name{PLIquantile_multivar}
\alias{PLIquantile_multivar}

\title{
  Perturbed-Law based sensitivity Indices (PLI) for quantile 
  and simultaneous perturbations of 2 inputs}

\description{
  PLIquantile_multivar computes the Perturbed-Law based Indices (PLI) for quantile
  and simultaneous perturbations of the means of 2 inputs, 
  estimated by a Monte Carlo method.
}

\usage{
PLIquantile_multivar(order,x,y,inputs,deltasvector,InputDistributions,samedelta=TRUE,
            percentage=TRUE,nboot=0,conf=0.95,bootsample=TRUE)
}

\arguments{
  \item{order}{the order of the quantile to estimate.}
  \item{x}{the matrix of simulation points coordinates, 
    one column per variable.}
    \item{y}{the vector of model outputs.}
  \item{inputs}{the vector of the two inputs' indices for which
    the indices will be computed.}  
	\item{deltasvector}{a vector containing the values of the perturbed means
	  for which the indices will be computed. Warning: if samedelta=FALSE, 
	  deltasvector has to be the vector of deltas (mean perturbations)}
	\item{InputDistributions}{a list of list. Each list contains, as a list,
    the name of the distribution to be used and the parameters.
	  Implemented cases so far (for a mean perturbation): 
	  Gaussian, Uniform, Triangle, Left Trucated Gaussian, 
	  Left Truncated Gumbel. Using Gumbel requires the package \code{evd}.}
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
      \item If it is set at FALSE, the initially proposed formula is used
      (see Sueur et al., 2017).
      \item If not (set as TRUE), the PLI is given in percentage of variation 
      of the quantile (see Iooss et al., 2021).
        }}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level for bootstrap confidence intervals.}
  \item{bootsample}{If TRUE, the uncertainty about the original quantile estimation
      is taken into account in the PLI confidence intervals (see Iooss et al., 2021). 
      If FALSE, standard confidence intervals are computed for the PLI.
      It mainly changes the CI at small delta values.}
}

\value{
  \code{PLIquantile_multivar} returns a list of matrices 
  (delta twist of input 1 (in rows) vs. delta twist of input 2 (in columns)) 
  containing the following components:
  \item{PLI}{the PLI.}
  \item{PLICIinf}{the bootstrap lower confidence interval values of the PLI.}
  \item{PLICIsup}{the bootstrap upper confidence interval values of the PLI.}
  \item{quantile}{the perturbed quantile.}
  \item{quantileCIinf}{the bootstrap lower confidence interval values of the perturbed quantile.}
  \item{quantileCIsup}{the bootstrap upper confidence interval values of the perturbed quantile.}
}

\details{
This function does not allow perturbations on the variance of the inputs' distributions.
}

\references{
T. Delage, R. Sueur and B. Iooss, 2018, \emph{Robustness analysis of epistemic uncertainties 
propagation studies in LOCA assessment thermal-hydraulic model}, 
ANS Best Estimate Plus Uncertainty International Conference (BEPU 2018), Lucca, Italy, May 13-19, 2018.

B. Iooss, V. Verges and V. Larget, 2022, \emph{BEPU robustness analysis via perturbed law-based 
sensitivity indices}, Proceedings of the Institution of Mechanical Engineers, 
Part O: Journal of Risk and Reliability, 236:855-865.

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
Bertrand Iooss
}

\seealso{
  \code{\link{PLI}, \link{PLIquantile}, \link{PLIsuperquantile}, \link{PLIsuperquantile_multivar}}
}

\examples{
\donttest{

# Model: 3D function 

distribution = list()
for (i in 1:3) distribution[[i]]=list("norm",c(0,1))
N = 5000
X = matrix(0,ncol=3,nrow=N)
for(i in 1:3) X[,i] = rnorm(N,0,1)
Y = 2 * X[,1] + X[,2] + X[,3]/2
alpha <- 0.95
nboot <- 20 # put nboot=200 for consistency

q95 = quantile(Y,alpha)
v_delta = seq(-1,1,1/10) 
toto12 = PLIquantile_multivar(alpha,X,Y,c(1,2),deltasvector=v_delta,
    InputDistributions=distribution,samedelta=TRUE)
toto = PLIquantile(alpha,X,Y,deltasvector=v_delta,InputDistributions=distribution,
    type="MOY",samedelta=TRUE,nboot=0)

par(mar=c(4,5,1,1))
plot(v_delta,diag(toto12$PLI),,ylim=c(-1,1),xlab=expression(delta),
    ylab=expression(hat(PLI[i*delta])),pch=16,cex=1.5,col="blue")
points(v_delta,toto$PLI[,1],col="darkgreen",pch=15,cex=1.5)
points(v_delta,toto$PLI[,2],col="black",pch=19,cex=1.5)
points(v_delta,toto$PLI[,3],col="red",pch=17,cex=1.5)
abline(h=0,lty=2)
legend(-1,1.,legend=c("X1","X2","X3","X1X2"),col=c("darkgreen","black","red","blue"),
    pch=c(15,19,17,16),cex=1.5)

# with bootstrap

v_delta = seq(-1,1,2/10) 

toto12 = PLIquantile_multivar(alpha,X,Y,c(1,2),deltasvector=v_delta,
    InputDistributions=distribution,samedelta=TRUE,nboot=nboot,bootsample=FALSE)
toto = PLIquantile(alpha,X,Y,deltasvector=v_delta,InputDistributions=distribution,
    type="MOY",samedelta=TRUE,nboot=nboot,bootsample=FALSE)

par(mar=c(4,5,1,1))
plot(v_delta,diag(toto12$PLI),ylim=c(-1,1),xlab=expression(delta),
    ylab=expression(hat(PLI[i*delta])),pch=16,cex=1.5,col="blue")
points(v_delta,toto$PLI[,1],col="darkgreen",pch=15,cex=1.5)
points(v_delta,toto$PLI[,2],col="black",pch=19,cex=1.5)
points(v_delta,toto$PLI[,3],col="red",pch=17,cex=1.5)
lines(v_delta,diag(toto12$PLICIinf),col="blue")
lines(v_delta,diag(toto12$PLICIsup),col="blue")
lines(v_delta,toto$PLICIinf[,2],col="black")
lines(v_delta,toto$PLICIsup[,2],col="black")
lines(v_delta,toto$PLICIinf[,1],col="darkgreen")
lines(v_delta,toto$PLICIsup[,1],col="darkgreen")
lines(v_delta,toto$PLICIinf[,3],col="red")
lines(v_delta,toto$PLICIsup[,3],col="red")
abline(h=0,lty=2)
legend(-1,1,legend=c("X1","X2","X3","X1X2"),col=c("darkgreen","black","red","blue"),
    pch=c(15,19,17,16),cex=1.5)

###################################################		
# another visualizations by using the plotrix, 
# viridisLite, lattice and grid packages (from Vanessa Verges)

library(plotrix)

parameters = list(colors=c("darkgreen","black","red"),symbols=c(15,19,17))
par(mar=c(4,5,1,1),xpd=TRUE)
plotCI(v_delta,diag(toto12$PLI),ui=diag(toto12$PLICIsup),li=diag(toto12$PLICIinf),
       xlab=expression(delta),ylab=expression(hat(PLI[i*delta])),
       main=bquote("PLI-quantile (N ="~.(N) ~ ","~alpha~"="~.(alpha)~
       ") on "~X[1]~"and"~X[2]~"of Y="~2*X[1] + X[2] + X[3]/2),
       cex=1.5,col="blue",pch=16)
for (i in 1:3){
  plotCI(v_delta,toto$PLI[,i],ui=toto$PLICIsup[,i],li=toto$PLICIinf[,i],
         cex=1.5,col=parameters$colors[i],pch=parameters$symbols[i],
         add=TRUE)
}
abline(h=0,lty=2)
legend("topleft",legend=c("X1","X2","X3","X1X2"),
        col=c(parameters$colors,"blue"),pch=c(parameters$symbols,16),cex=1.5)

# Visu of all the PLIs (at any paired combinations of deltas)

library(viridisLite)
library(lattice)
library(grid)

colnames(toto12$PLI) = round(v_delta,2)
rownames(toto12$PLI) = round(v_delta,2)
coul = viridis(100)
levelplot(toto12$PLI, col.regions = coul, xlab=bquote(delta[X~.(1)]), ylab=bquote(delta[X~.(2)]),
         main=bquote(hat(PLI)[quantile[~X[1]~X[2]]]))

}
}
