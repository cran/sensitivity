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
    \item Linear regression coefficients: SRC and SRRC (\code{\link{src}}), PCC and PRCC (\code{\link{pcc}});
    \item Bettonvil's sequential bifurcations (Bettonvil and Kleijnen, 1996) (\code{\link{sb}});
    \item Morris's "OAT" elementary effects screening method (\code{\link{morris}});
    \item Derivative-based Global Sensitivity Measures:
    \itemize{
      \item Poincare constants for Derivative-based Global Sensitivity Measures (DGSM) (Lambon et al., 2013; Roustant et al., 2016) (\code{\link{PoincareConstant}}) and (\code{\link{PoincareOptimal}}),
      \item Distributed Evaluation of Local Sensitivity Analysis (DELSA) (Rakovec et al., 2014) (\code{\link{delsa}});
      }
    \item Variance-based sensitivity indices (Sobol' indices):
	  \itemize{
      \item Estimation of the Sobol' first order indices with with B-spline Smoothing (Ratto and Pagano, 2010) (\code{\link{sobolSmthSpl}}),
      \item Monte Carlo estimation of Sobol' indices with independent inputs (also called pick-freeze method): 
      \itemize{
		    \item Sobol' scheme (Sobol, 1993) to compute the indices given by the variance decomposition up to a specified order (\code{\link{sobol}}),
  		  \item Saltelli's scheme (Saltelli, 2002) to compute first order, second order and total indices  (\code{\link{sobolSalt}}), 
  		  \item Saltelli's scheme (Saltelli, 2002) to compute first order and total indices (\code{\link{sobol2002}}), 
	  	  \item Mauntz-Kucherenko's scheme (Sobol et al., 2007) to compute first order and total indices using improved formulas for small indices (\code{\link{sobol2007}}),
	  	  \item Jansen-Sobol's scheme (Jansen, 1999) to compute first order and total indices using improved formulas (\code{\link{soboljansen}}),
    	  \item Martinez's scheme using correlation coefficient-based formulas (Martinez, 2011; touati, 2016) to compute first order and total indices, associated with theoretical confidence intervals (\code{\link{sobolmartinez}} and \code{\link{soboltouati}}), 
	  	  \item Janon-Monod's scheme (Monod et al., 2006; Janon et al., 2013) to compute first order indices with optimal asymptotic variance (\code{\link{sobolEff}}),
	  	  \item Mara's scheme (Mara and Joseph, 2008) to compute first order indices with a cost independent of the dimension, via a unique-matrix permutations (\code{\link{sobolmara}}),
		    \item Owen's scheme (Owen, 2013) to compute first order and total indices using improved formulas (via 3 input independent matrices) for small indices (\code{\link{sobolowen}}),
        \item Total Interaction Indices using Liu-Owen's scheme (Liu and Owen, 2006) (\code{\link{sobolTIIlo}}) and pick-freeze scheme (Fruth et al., 2014) (\code{\link{sobolTIIpf}}),
      }
      \item Estimation of the Sobol' first order and total indices with Saltelli's so-called "extended-FAST" method (Saltelli et al., 1999) (\code{\link{fast99}}),
      \item Estimation of the Sobol' first order and closed second order indices using replicated orthogonal array-based Latin hypecube sample (Tissot and Prieur, 2015) (\code{\link{sobolroalhs}}),
      \item Sobol' indices estimation under inequality constraints (Gilquin et al., 2015) by extension of the replication procedure (Tissot and Prieur, 2015) (\code{\link{sobolroauc}}),
      \item Estimation of the Sobol' first order and total indices with kriging-based global sensitivity analysis (Le Gratiet et al., 2014) (\code{\link{sobolGP}});
    }
    \item Variance-based sensitivity indices (Shapley effects and Sobol' indices, with independent or dependent inputs):
	  \itemize{
	    \item Estimation by examining all permutations of inputs (Song et al., 2016) (\code{\link{shapleyPermEx}})
	    \item Estimation by randomly sampling permutations of inputs (Song et al., 2016) (\code{\link{shapleyPermRand}})
	  }
    \item Support index functions (\code{\link{support}}) of Fruth et al. (2015);
    \item Sensitivity Indices based on Csiszar f-divergence (\code{\link{sensiFdiv}}) (particular cases: Borgonovo's indices and mutual-information based indices) and Hilbert-Schmidt Independence Criterion (\code{\link{sensiHSIC}}) of Da Veiga (2015);
    \item Reliability sensitivity analysis by the Perturbed-Law based Indices (\code{\link{PLI}}) of Lemaitre et al. (2015) and (\code{\link{PLIquantile}}) of Sueur et al. (2016);
    \item Sobol' indices for multidimensional outputs (\code{\link{sobolMultOut}}): Aggregated Sobol' indices (Lamboni et al., 2011; Gamboa et al., 2014) and functional (1D) Sobol' indices.
  }
  
  Moreover, some utilities are provided: standard test-cases
  (\code{\link{testmodels}}) and template file generation
  (\code{\link{template.replace}}).
}

\section{Model managing}{
  The \pkg{sensitivity} package has been designed to work either models written in \R
  than external models such as heavy computational codes. This is achieved with
  the input argument \code{model} present in all functions of this package.
  
  The argument \code{model} is expected to be either a
  funtion or a predictor (i.e. an object with a \code{predict} function such as
  \code{lm}).
  \itemize{
    \item If \code{model = m} where \code{m} is a function, it will be invoked
        once by \code{y <- m(X)}.
    \item If \code{model = m} where \code{m} is a predictor, it will be invoked
        once by \code{y <- predict(m, X)}.
  }
  \code{X} is the design of experiments, i.e. a \code{data.frame} with
  \code{p} columns (the input factors) and \code{n} lines (each, an
  experiment), and \code{y} is the vector of length \code{n} of the
  model responses.
  
  The model in invoked once for the whole design of experiment.
  
  The argument \code{model} can be left to \code{NULL}. This is refered to as 
  the decoupled approach and used with external computational codes that rarely
  run on the statistician's computer. See \code{\link{decoupling}}.
}


\author{Gilles Pujol, Bertrand Iooss, Alexandre Janon with contributions from Paul Lemaitre for the \code{\link{PLI}} function, Laurent Gilquin for the \code{\link{sobolroalhs}}, \code{\link{sobolroauc}} and \code{\link{sobolSalt}} functions, Loic le Gratiet for the \code{\link{sobolGP}} function, Khalid Boumhaout, Taieb Touati and Bernardo Ramos for the \code{\link{sobolowen}} and \code{\link{soboltouati}} functions, Jana Fruth for the \code{\link{PoincareConstant}}, \code{\link{support}}, \code{\link{sobolTIIlo}} and \code{\link{sobolTIIpf}} functions, Sebastien Da veiga for the \code{\link{sensiFdiv}} and \code{\link{sensiHSIC}} functions, Joseph Guillaume for the \code{\link{delsa}} and \code{\link{parameterSets}} functions, Olivier Roustant for the \code{\link{PoincareOptimal}} function, Eunhye Song, Barry L. Nelson and Jeremy Staum for the \code{\link{shapleyPermEx}} and \code{\link{shapleyPermRand}} functions, Filippo Monari for the (\code{\link{sobolSmthSpl}}) function, Frank Weber, Thibault Delage and Roelof Oomen.

(maintainer: Bertrand Iooss \email{biooss@yahoo.fr})}

\references{
R. Faivre, B. Iooss, S. Mahevas, D. Makowski, H. Monod, editors, 2013, \emph{Analyse de sensibilite et exploration de modeles. Applications aux modeles environnementaux}, Editions Quae.

B. Iooss and A. Saltelli, 2017, \emph{Introduction: Sensitivity analysis.} In: \emph{Springer Handbook on Uncertainty Quantification}, R. Ghanem, D. Higdon and H. Owhadi (Eds), Springer. href{http://link.springer.com/referenceworkentry/10.1007/978-3-319-11259-6_31-1}

B. Iooss and P. Lemaitre, 2015, \emph{A review on global sensitivity analysis methods}. In \emph{Uncertainty management in Simulation-Optimization of Complex Systems: Algorithms and Applications}, C. Meloni and G. Dellino (eds), Springer.

A. Saltelli, K. Chan and E. M. Scott eds, 2000, \emph{Sensitivity Analysis}, Wiley.

A. Saltelli et al., 2008, \emph{Global Sensitivity Analysis: The Primer}, Wiley

}

\keyword{package}

