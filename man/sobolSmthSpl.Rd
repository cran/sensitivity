%Documentation for sobolSmthSpl.
\name{sobolSmthSpl}
\alias{sobolSmthSpl}
\title{Estimation of Sobol' First Order Indices with B-spline Smoothing}
\description{
	Determines the Si coefficient for singular parameters through B-spline smoothing with roughness penalty.
}
\usage{
	sobolSmthSpl(Y, X)
}
\arguments{
	\item{Y}{vector of model responses.}
	\item{X}{matrix having as rows the input vectors corresponding to the responses in Y.}	
}
\value{
	sobolSmthSpl returns a list of class "sobolSmthSpl" containing the following components:
	\item{call}{the matched call.}
	\item{X}{the provided input matrix.}
	\item{Y}{the provided matrix of model responses.}
	\item{S}{a matrix having the following columns:
		Si (the estimated first order Sobol' indices), 
		Si.e (the standard errors for the estimated first order Sobol' indices) and
		q0.05 (the 0.05 quantiles assuming for the Si indices Normal distributions centred on the 
			Si estimates and with standard deviations the calculated standard errors)
	}
}
\details{
WARNING: This function can give bad results for reasons that have not been yet investigated.
}
\author{
	Filippo Monari
}
\references{
	Saltelli, A; Ratto, M; Andres, T; Campolongo, F; Cariboni, J; Gatelli, D; Saisana, M & Tarantola, S.
	\emph{Global Sensitivity Analysis: The Primer Wiley-Interscience}, 2008
	
  M Ratto and A. Pagano, 2010, \emph{Using recursive algorithms for the efficient identification
  of smoothing spline ANOVA models}, Advances in Statistical Analysis, 94, 367--388.
	
}
\seealso{
	 \code{\link{sobol}, \link{sobolEff}, \link{sobolGP}}
}

\examples{
	X = matrix(runif(5000), ncol = 10)
	Y = sobol.fun(X)
	sa = sobolSmthSpl(Y, X)
	plot(sa)
}
