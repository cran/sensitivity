\name{sobolCert}
\alias{sobolCert}
\alias{print.sobolCert}

\title{Monte Carlo Estimation of Sobol' Indices using certified meta-models}

\description{ \code{sobolCert} implements the Monte Carlo estimation of
   the Sobol' sensitivity indices using certified metamodels using
   the formulas in Janon et al. (2011). }

\usage{
sobolCert(model = NULL, X1=NULL, X2=NULL, nboot = 300, conf = 0.95, lambda0 = 0, h = 0)
\method{print}{sobolCert}(x, \dots)
}

\arguments{
  \item{model}{a function defining the model to analyze. This function must return a list whose components are:
  \itemize{
	  \item{out}{metamodel output.}
	  \item{err}{metamodel output error bound, satisfying \deqn{|model_{output} - metamodel_{output}| <= err}}
	  }}
  \item{X1}{the first random sample. If NULL, sobolEff ignores model, X1 and X2, and reuse model outputs from the previous sobolCert call. }
  \item{X2}{the second random sample.}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level for confidence intervals.}
  \item{lambda0}{if lambda0=0, use the method described in Janon et al. (2011) Section 3.1; else use the method of Section 3.2 with lambda0 as penalty parameter.}
  \item{h}{if lambda0=0, this parameter is ignored; else it is used as the h bandwidth parameter.}
  \item{x}{a sobolCert object}
  \item{...}{currently not used}
}

\value{
  \code{sobolCert} returns a list of class \code{"sobolCert"}, containing the following components:

  \item{call}{the matched call.}
  \item{S}{the estimations of the Sobol' sensitivity indices.}
  \item{penalty}{(only if lambda0>0) value of the smoothing penalty.}
}

\references{
	Janon, A., Nodet M., Prieur C. (2011) Uncertainties assessment in global 
  sensitivity indices estimation from metamodels. 
  To appear in International Journal for Uncertainty Quantification. 
}

\author{
Alexandre Janon
}

\seealso{
  \code{\link{sobol}, \link{sobol2002}, \link{sobol2007}}
}

\examples{
	\dontrun{
# Test case

n <- 1000
X1 <- data.frame(matrix(runif(3 * n), nrow = n))
X2 <- data.frame(matrix(runif(3 * n), nrow = n))

# sensitivity analysis
x=sobolCert(model=function(X) { list(out=X[1]+2*X[2]+X[3]+.001*runif(1),err=.01); }, 
            X1, X2, conf=.99, lambda0=.1, h=.1, nboot=30)
print(x)

x=sobolCert(model=NULL, X1=NULL, X2=NULL, conf=.95)
print(x)
}
}

\keyword{design}
