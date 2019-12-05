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
            percentage=FALSE)
}

\arguments{
  \item{order}{the order of the quantile to estimate.}
  \item{x}{the matrix of simulation points coordinates, 
    one column per variable.}
    \item{y}{the vector of model outputs.}
  \item{inputs}{the vector of the two inputs' indices for which
    the indices will be computed.}  
	\item{deltasvector}{a vector containing the values of delta for which
    the indices will be computed.}
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
      \item If it is set at FALSE, the classical formula used in the bibliographical 
      references is used.
      \item If not (set as TRUE), the PLI is given in percentage of variation 
      of the quantile (even if it is negative).
        }}
}

\value{
  \code{PLIquantile_multivar} returns a list of matrix (each column corresponds to an input, 
  each line corresponds to a twist of amplitude delta) 
  containing the following components:
  \item{PLI}{the PLI.}
  \item{quantile}{the perturbed quantile.}
}

\details{
This function does not allow perturbations on the variance of the inputs' distributions.
This function does not allow bootstrap in order to obtain confidence intervals on the PLI estimates.
}

\references{
T. Delage, R. Sueur and B. Iooss, 2018, \emph{Robustness analysis of epistemic uncertainties 
propagation studies in LOCA assessment thermal-hydraulic model}, 
ANS Best Estimate Plus Uncertainty International Conference (BEPU 2018), Lucca, Italy, May 13-19, 2018.

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
  \code{\link{PLI}, \link{PLIquantile}, \link{PLIsuperquantile}}
}

\examples{
\donttest{

# Model: 3D function 

distribution = list()
for (i in 1:3) distribution[[i]]=list("norm",c(0,1))
N = 10000
X = matrix(0,ncol=3,nrow=N)
for(i in 1:3) X[,i] = rnorm(N,0,1)
Y = 2 * X[,1] + X[,2] + X[,3]/2

q95 = quantile(Y,0.95)
v_delta = seq(-1,1,1/10) 
toto12 = PLIquantile_multivar(0.95,X,Y,c(1,2),deltasvector=v_delta,
    InputDistributions=distribution,samedelta=TRUE)
toto = PLIquantile(0.95,X,Y,deltasvector=v_delta,InputDistributions=distribution,
    type="MOY",samedelta=TRUE,nboot=0)

par(mar=c(4,5,1,1))
plot(v_delta,diag(toto12$PLI),ylim=c(-2.5,1.5),xlab=expression(delta),
    ylab=expression(hat(PLI[i*delta])),pch=16,cex=1.5,col="blue")
points(v_delta,toto$PLI[,1],col="darkgreen",pch=15,cex=1.5)
points(v_delta,toto$PLI[,2],col="black",pch=19,cex=1.5)
points(v_delta,toto$PLI[,3],col="red",pch=17,cex=1.5)
abline(h=0,lty=2)
legend(-1,1.5,legend=c("X1","X2","X3","X1X2"),col=c("darkgreen","black","red","blue"),
    pch=c(15,19,17,16),cex=1.5)

}
}
