\name{sobolGP}

\alias{sobolGP}
\alias{plot.sobolGP}
\alias{ask.sobolGP}
\alias{tell.sobolGP}
\alias{print.sobolGP}

\title{Kriging-based sensitivity analysis}

\description{
Perform a kriging-based global sensitivity analysis taking into account both the meta-model and the Monte-Carlo errors. The Sobol indices are estimated with a Monte-Carlo integration and the true function is substituted by a kriging model. It is built thanks to the function \code{km} of the package \code{DiceKriging}.
The complete conditional predictive distribution of the kriging model is considered (not only the predictive mean).
}
\usage{
sobolGP(
	model,
	type="SK",
	MCmethod="sobol",
	X1,
	X2, 
	nsim=100,
	nboot=1,
	conf = 0.95,
	sequential = FALSE, 
	candidate, 
	sequential.tot=FALSE,
	max_iter = 1000) 

\method{ask}{sobolGP}(x, tot = FALSE, \dots)

\method{tell}{sobolGP}(x,y=NULL,xpoint=NULL,newcandidate=NULL, \dots)

\method{print}{sobolGP}(x, \dots)

\method{plot}{sobolGP}(x,\dots)
}

\arguments{
  \item{model}{
an object of class \code{"km"} specifying the kriging model built from package \code{"DiceKriging"} (see \code{\link{km}}).
}
  \item{type}{
a character string giving the type of the considered kriging model. \code{"SK"} refers to simple kriging and \code{"UK"} refers to universal kriging (see \code{\link{km}}).
}
  \item{MCmethod}{
a character string specifying the Monte-Carlo procedure used to estimate the Sobol indices. The avaible methods are : \code{"sobol"}, \code{"sobol2002"}, \code{"sobol2007"}, \code{"sobolEff"} and \code{"soboljansen"}.
}
  \item{X1}{
a matrix representing the first random sample.  
}
  \item{X2}{
a matrix representing the second random sample.
}
  \item{nsim}{
an integer giving the number of samples for the conditional Gaussian process. It is used to quantify the uncertainty due to the kriging approximation.
}
  \item{nboot}{
an integer representing the number of bootstrap replicates. It is used to quantify the uncertainty due to the Monte-Carlo integrations. We recommend to set \code{nboot = 100}.
}
  \item{conf}{
a numeric representing the confidence intervals taking into account the uncertainty due to the bootstrap procedure and the Gaussian process samples.
}
  \item{sequential}{
a boolean. If \code{sequential=TRUE}, the procedure provides a new point where to perform a simulation. It is the one minimizing the sum of the MAIN effect estimate variances. The variance is taken with respect to the conditional Gaussian process. The new point is selected in the points \code{candidate}.
}
  \item{candidate}{
a matrix representing the candidate points where the best new point to be simulated is selected. The lines represent the points and the columns represent the dimension.
}
  \item{sequential.tot}{
a boolean. If \code{sequential.tot=TRUE}, the procedure provides a new point where to perform the simulation. It is the one minimizing the sum of the TOTAL effect estimate. The variance is taken with respect to the conditional Gaussian process. The new point is selected in the points \code{candidate}.
}
  \item{max_iter}{
a numeric giving the maximal number of iterations for the propagative Gibbs sampler. It is used to simulate  the realizations of the Gaussian process.
}
  \item{x}{
an object of class S3 \code{"sobolGP"} obtaining from the procedure \code{sobolGP}. It stores the results of the Kriging-based global sensitivity analysis.
}
  \item{tot}{
a boolean. If \code{tot=TRUE}, the procedure ask provides a point relative to the uncertainty of the total Sobol' indices (instead of first order' ones). 
}
  \item{xpoint}{
a matrix representing a new point added to the kriging model.
}
  \item{y}{
a numeric giving the response of the function at \code{xpoint}.
}
  \item{newcandidate}{
a matrix representing the new candidate points where the best point to be simulated is selected. If \code{newcandidate=NULL}, these points correspond to \code{candidate} without the new point \code{xpoint}.
  }
  \item{\dots}{any other arguments to be passed}
}
\value{
An object of class S3 \code{sobolGP}.
\itemize{
\item call : a list containing the arguments of the function \code{sobolGP} :
	\itemize{
	\item X1 : X1
	\item X2 : X2
	\item conf : conf
	\item nboot : nboot
	\item candidate : candidate
	\item sequential : sequential
	\item max_iter : max_iter
	\item sequential.tot : sequential.tot
	\item model : model
	\item tot : tot
	\item method : MCmethod
	\item type : type
	\item nsim : nsim
	}
\item S : a list containing the results of the kriging-based sensitivity analysis for the MAIN effects:
	\itemize{
	\item mean : a matrix giving the mean of the Sobol index estimates.
	\item var : a matrix giving the variance of the Sobol index estimates.
	\item ci : a matrix giving the confidence intervals of the Sobol index estimates according to \code{conf}.
	\item varPG : a matrix giving the variance of the Sobol index estimates due to the Gaussian process approximation.
	\item varMC : a matrix giving the variance of the Sobol index estimates due to the Monte-Carlo integrations.
	\item xnew : if \code{sequential=TRUE}, a matrix giving the point in \code{candidate} which is the best to simulate.
	\item xnewi : if \code{sequential=TRUE}, an integer giving the index of the point in \code{candidate} which is the best to simulate.
	}
\item T : a list containing the results of the kriging-based sensitivity analysis for the TOTAL effects:
	\itemize{
	\item mean : a matrix giving the mean of the Sobol index estimates.
	\item var : a matrix giving the variance of the Sobol index estimates.
	\item ci : a matrix giving the confidence intervals of the Sobol index estimates according to \code{conf}.
	\item varPG : a matrix giving the variance of the Sobol index estimates due to the Gaussian process approximation.
	\item varMC : a matrix giving the variance of the Sobol index estimates due to the Monte-Carlo integrations.
	\item xnew : if \code{sequential.tot=TRUE}, a matrix giving the point in \code{candidate} which is the best to simulate.
	\item xnewi : if \code{sequential.tot=TRUE}, an integer giving the index of the point in \code{candidate} which is the best to simulate.
	}
}
}

\details{
	The function \code{ask} provides the new point where the function should be simulated. Furthermore, the function \code{tell} performs a new kriging-based sensitivity analysis when the point \code{x} with the corresponding observation \code{y} is added.
}

\references{
Le gratiet, L., Cannamela, C. and Iooss, B. (2014), A Bayesian approach for global sensitivity analysis of (multi-fidelity) computer codes, \emph{ arXiv:1307.2223}
}

\author{
Loic Le Gratiet, EDF R&D - CNRS, I3S
}


\seealso{
\code{\link{sobol}}, \code{\link{sobol2002}}, \code{\link{sobol2007}}, \code{\link{sobolEff}}, \code{\link{soboljansen}}, \code{\link{km}}
}
\examples{

\dontrun{
library(DiceKriging)

#--------------------------------------#
# kriging model building
#--------------------------------------#

d <- 2; n <- 16
design.fact <- expand.grid(x1=seq(0,1,length=4), x2=seq(0,1,length=4))
y <- apply(design.fact, 1, branin) 

m <- km(design=design.fact, response=y)

#--------------------------------------#
# sobol samples & candidate points
#--------------------------------------#

n <- 1000
X1 <- data.frame(matrix(runif(d * n), nrow = n))
X2 <- data.frame(matrix(runif(d * n), nrow = n))

candidate <- data.frame(matrix(runif(d * 100), nrow = 100))

#--------------------------------------#
# Kriging-based Sobol
#--------------------------------------#

res <- sobolGP(
model = m,
type="UK",
MCmethod="sobol",
X1,
X2,
nsim = 100,
conf = 0.95,
nboot=100,
sequential = TRUE,
candidate,
sequential.tot=FALSE,
max_iter = 1000
) 

res
plot(res)
x <- ask.sobolGP(res)
y <- branin(x)
res.new <- tell.sobolGP(res,y,x)
res.new
}
}

