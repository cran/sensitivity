\name{sensitivity-package}
\alias{sensitivity-package}
\alias{sensitivity}
\docType{package}

\title{Sensitivity Analysis}

\description{Methods and functions for global sensitivity analysis.}

\details{
  The \pkg{sensitivity} package implements some global sensitivity analysis
  methods:
  \itemize{
    \item Linear regression coefficients: SRC and SRRC
    (\code{\link{src}}), PCC and PRCC (\code{\link{pcc}}).
    \item Morris's "OAT" elementary effects screening method (\code{\link{morris}}).
    \item Bettonvil's sequential bifurcations (\code{\link{sb}}).
    \item Monte Carlo estimation of Sobol' indices: 
	  \itemize{
		  \item Sobol's scheme (1993) to compute the indices given by the variance decomposition up to a specified order (\code{\link{sobol}})
  		  \item Saltelli's scheme (2002) to compute first order and total indices with a reduced cost (\code{\link{sobol2002}}).
		  \item Monod's scheme (2006) to compute first order indices with optimal asymptotic variance (\code{\link{sobolEff}}). 
  	  \item Sobol's scheme (2007) to compute first order and total indices using improved formulas for small indices (\code{\link{sobol2007}}).
  	  \item Jansen-Sobol's scheme (1999) to compute first order and total indices using improved formulas (\code{\link{soboljansen}}).
	  }
    \item Estimation of the Sobol' first order and total indices with
    Saltelli's so-called "extended-FAST" method (\code{\link{fast99}}).
  }
  
  Moreover, some utilities are provided: standard test-cases
  (\code{\link{testmodels}}) and template file generation
  (\code{\link{template.replace}}).
}

\section{Model managing}{
  The \pkg{sensitivity} package has been designed to work either models written in \R
  than external models such as heavy computational codes. This is achieved with
  the input argument \code{model} present in all functions of this package.
  
  If \code{model} isn't set to \code{NULL}, it is expected to be either a
  funtion or a predictor (i.e. an object with a \code{predict} function such as
  \code{lm}).
  \itemize{
    \item If \code{model = f} where \code{f} is a function, it will be invoked
        once by \code{y <- f(X)}.
    \item If \code{model = r} where \code{r} is a predictor, it will be invoked
        once by \code{y <- predict(r, X)}.
  }
  \code{X} is the design of experiments, i.e. a \code{data.frame} with
  \code{p} columns (the input factors) and \code{n} lines (each, an
  experiment), and \code{y} is the vector of length \code{n} of the
  model responses.
  
  We draw attention to the fact that the model in invoked once for the whole
  design of experiment, which can be unnatural to many users. If the model to analyze
  isn't vectorized, it can be encapsulated in a parent loop model. It can also been used
  as it through the decoupled approach as described in \code{\link{decoupling}}.
  This is also how one can use the functions of the package with external
  computational codes that rarely run on the statistician's computer.
}

\references{
A. Saltelli, K. Chan and E. M. Scott eds, 2000, \emph{Sensitivity
  Analysis}, Wiley.
}

\keyword{package}

