\name{sensitivity-package}
\alias{sensitivity-package}
\alias{sensitivity}
\docType{package}

\title{Sensitivity Analysis}

\description{Methods and functions for global sensitivity analysis}

\details{
  The \pkg{sensitivity} package implements some global sensitivity analysis
  methods:
  \itemize{
    \item Linear regression importance measures (work in regression and classification contexts):
  \itemize{
    \item SRC and SRRC (\code{\link{src}}), 
    \item PCC, SPCC, PRCC and SPRCC (\code{\link{pcc}}),
    \item LMG and LMG on ranks (\code{\link{lmg}})
    \item E-MVD and E-MVD on ranks (\code{\link{lmg}});
    }
    \item Bettonvil's sequential bifurcations (Bettonvil and Kleijnen, 1996) (\code{\link{sb}});
    \item Morris's "OAT" elementary effects screening method (\code{\link{morris}});
    \item Derivative-based Global Sensitivity Measures:
    \itemize{
      \item Poincare constants for Derivative-based Global Sensitivity Measures (DGSM) (Lamboni et al., 2013; Roustant et al., 2017) (\code{\link{PoincareConstant}}) and (\code{\link{PoincareOptimal}}),
      \item Squared coefficients computation in generalized chaos via Poincare differential operators (Roustant et al., 2019) (\code{\link{PoincareChaosSqCoef}}),
      \item Distributed Evaluation of Local Sensitivity Analysis (DELSA) (Rakovec et al., 2014) (\code{\link{delsa}});
      }
    \item Variance-based sensitivity indices (Sobol' indices) for independent inputs:
	  \itemize{
      \item Estimation of the Sobol' first order indices with with B-spline Smoothing (Ratto and Pagano, 2010) (\code{\link{sobolSmthSpl}}),
      \item Monte Carlo estimation of Sobol' indices with independent inputs (also called pick-freeze method): 
      \itemize{
		    \item Sobol' scheme (Sobol, 1993) to compute the indices given by the variance decomposition up to a specified order (\code{\link{sobol}}),
  		  \item Saltelli's scheme (Saltelli, 2002) to compute first order, second order and total indices  (\code{\link{sobolSalt}}), 
  		  \item Saltelli's scheme (Saltelli, 2002) to compute first order and total indices (\code{\link{sobol2002}}), 
	  	  \item Mauntz-Kucherenko's scheme (Sobol et al., 2007) to compute first order and total indices using improved formulas for small indices (\code{\link{sobol2007}}),
	  	  \item Jansen-Sobol's scheme (Jansen, 1999) to compute first order and total indices using improved formulas (\code{\link{soboljansen}}),
    	  \item Martinez's scheme using correlation coefficient-based formulas (Martinez, 2011; Touati, 2016) to compute first order and total indices, associated with theoretical confidence intervals (\code{\link{sobolmartinez}} and \code{\link{soboltouati}}), 
	  	  \item Janon-Monod's scheme (Monod et al., 2006; Janon et al., 2013) to compute first order indices with optimal asymptotic variance (\code{\link{sobolEff}}),
	  	  \item Mara's scheme (Mara and Joseph, 2008) to compute first order indices with a cost independent of the dimension, via permutations on a single matrix (\code{\link{sobolmara}}),
	  	  \item Mighty estimator of first-order sensitivity indices based on rank statistics (correlation coefficient of Chatterjee, 2019; Gamboa et al., 2020) (\code{\link{sobolrank}}),
		    \item Owen's scheme (Owen, 2013) to compute first order and total indices using improved formulas (via 3 input independent matrices) for small indices (\code{\link{sobolowen}}),
        \item Total Interaction Indices using Liu-Owen's scheme (Liu and Owen, 2006) (\code{\link{sobolTIIlo}}) and pick-freeze scheme (Fruth et al., 2014) (\code{\link{sobolTIIpf}}),
      }
      \item Replication-based procedures:
      \itemize{
        \item Estimation of the Sobol' first order and closed second order indices using replicated orthogonal array-based Latin hypecube sample (Tissot and Prieur, 2015) (\code{\link{sobolroalhs}}),
        \item Recursive estimation of the Sobol' first order and closed second order indices using replicated orthogonal array-based Latin hypecube sample (Gilquin et al., 2016) (\code{\link{sobolrec}}),
        \item Estimation of the Sobol' first order, second order and total indices using the generalized method with replicated orthogonal array-based Latin hypecube sample (Tissot and Prieur, 2015) (\code{\link{sobolrep}}),
        \item Sobol' indices estimation under inequality constraints (Gilquin et al., 2015) by extension of the replication procedure (Tissot and Prieur, 2015) (\code{\link{sobolroauc}}),
      }
      \item Estimation of the Sobol' first order and total indices with Saltelli's so-called "extended-FAST" method (Saltelli et al., 1999) (\code{\link{fast99}}),
      \item Estimation of the Sobol' first order and total indices with kriging-based global sensitivity analysis (Le Gratiet et al., 2014) (\code{\link{sobolGP}});
    }
    \item Variance-based sensitivity indices (Shapley effects and Sobol' indices) for independent or dependent inputs:
	  \itemize{
	    \item Exact computation in the linear Gaussian framework (Broto et al., 2019) (\code{\link{shapleyLinearGaussian}}),
	    \item Computation in the Gaussian linear framework with an unknown block-diagonal covariance matrix (Broto et al., 2020) (\code{\link{shapleyBlockEstimation}}),
	    \item Estimation by examining all permutations of inputs (Song et al., 2016) (\code{\link{shapleyPermEx}}),
	    \item Estimation by randomly sampling permutations of inputs (Song et al., 2016) (\code{\link{shapleyPermRand}}),
	    \item Estimation of Shapley effects from data using nearest neighbors method (Broto et al., 2018) (\code{\link{shapleySubsetMc}}),
	    \item Estimation of Shapley effects and all Sobol indices from data using nearest neighbors (Broto et al., 2018) using a fast approximate algorithm, and ranking (Gamboa et al., 2020) (\code{\link{shapleysobol_knn}}) and (\code{\link{sobolshap_knn}}),
	    \item Estimation of Shapley effects from data using nearest neighbors method (Broto et al., 2018) with an optimized/parallelized computations and bootstrap confidence intervals estimations  (\code{\link{shapleysobol_knn}});
	  }
    \item Support index functions (\code{\link{support}}) of Fruth et al. (2016);
    \item Sensitivity Indices based on Csiszar f-divergence (\code{\link{sensiFdiv}}) (particular cases: Borgonovo's indices and mutual-information based indices) and Hilbert-Schmidt Independence Criterion (\code{\link{sensiHSIC}}) (Da Veiga, 2015; Meynaoui et al., 2019);
    \item First-order quantile-oriented sensitivity indices as defined in Fort et al. (2016) via a kernel-based estimator related (Maume-Deschamps and Niang, 2018) (\code{\link{qosa}});
    \item Target Sensitivity Analysis via Hilbert-Schmidt Independence Criterion (\code{\link{sensiHSIC}}) (Spagnol et al., 2019);
    \item Robustness analysis by the Perturbed-Law based Indices (\code{\link{PLI}}) of Lemaitre et al. (2015), (\code{\link{PLIquantile}}) of Sueur et al. (2017), (\code{\link{PLIsuperquantile}}) of Iooss et al. (2021), and extension as (\code{\link{PLIquantile_multivar}}) and (\code{\link{PLIsuperquantile_multivar}}) ;
    \item Extensions to multidimensional outputs for:
	  \itemize{
	    \item Sobol' indices (\code{\link{sobolMultOut}}): Aggregated Sobol' indices (Lamboni et al., 2011; Gamboa et al., 2014) and functional (1D) Sobol' indices,
	    \item Shapley effects and Sobol' indices (\code{\link{shapleysobol_knn}}) and (\code{\link{sobolshap_knn}}): Functional (1D) indices,
	    \item HSIC indices (\code{\link{sensiHSIC}}) (Da Veiga, 2015): Aggregated HSIC, potentially via a PCA step (Da Veiga, 2015),
	    \item Morris method (\code{\link{morrisMultOut}}).
	    }
  }
  
  Moreover, some utilities are provided: standard test-cases (\code{\link{testmodels}}), weight transformation function of the output sample (\code{\link{weightTSA}}) to perform Target Sensitivity Analysis, normal and Gumbel truncated distributions (\code{\link{truncateddistrib}}), squared integral estimate (\code{\link{squaredIntEstim}}), Addelman and Kempthorne construction of orthogonal arrays of strength two (\code{\link{addelman_const}}), discrepancy criteria (\code{\link{discrepancyCriteria_cplus}}), maximin criteria (\code{\link{maximin_cplus}}) and template file generation (\code{\link{template.replace}}).
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


\author{Bertrand Iooss, Sebastien Da Veiga, Alexandre Janon and Gilles Pujol with contributions from Paul Lemaitre for \code{\link{PLI}}, Thibault Delage and Roman Sueur for \code{\link{PLIquantile}}, Vanessa Verges for \code{\link{PLIquantile}}, \code{\link{PLIsuperquantile}}, \code{\link{PLIquantile_multivar}} and \code{\link{PLIsuperquantile_multivar}}, Laurent Gilquin for \code{\link{sobolroalhs}}, \code{\link{sobolroauc}}, \code{\link{sobolSalt}}, \code{\link{sobolrep}}, \code{\link{sobolrec}}, as well as \code{\link{addelman_const}}, \code{\link{discrepancyCriteria_cplus}} and \code{\link{maximin_cplus}}, Loic le Gratiet for \code{\link{sobolGP}}, Khalid Boumhaout, Taieb Touati and Bernardo Ramos for \code{\link{sobolowen}} and \code{\link{soboltouati}}, Jana Fruth for \code{\link{PoincareConstant}}, \code{\link{sobolTIIlo}} and \code{\link{sobolTIIpf}}, Amandine Marrel, Anouar Meynaoui and Reda El Amri for their contributions to \code{\link{sensiHSIC}}, Joseph Guillaume and Oldrich Rakovec for \code{\link{delsa}} and \code{\link{parameterSets}}, Olivier Roustant for \code{\link{PoincareOptimal}}, \code{\link{PoincareChaosSqCoef}}, \code{\link{squaredIntEstim}} and \code{\link{support}}, Eunhye Song, Barry L. Nelson and Jeremy Staum for \code{\link{shapleyPermEx}} and \code{\link{shapleyPermRand}}, Baptiste Broto for \code{\link{shapleySubsetMc}}, \code{\link{shapleyLinearGaussian}} and \code{\link{shapleyBlockEstimation}}, Filippo Monari for (\code{\link{sobolSmthSpl}}) and (\code{\link{morrisMultOut}}), Marouane Il Idrissi for \code{\link{lmg}}, \code{\link{emvd}} and \code{\link{shapleysobol_knn}}, Frank Weber and Roelof Oomen for other contributions.

(maintainer: Bertrand Iooss \email{biooss@yahoo.fr})}

\references{
S. Da Veiga, F. Gamboa, B. Iooss and C. Prieur, \emph{Basics and trends in sensitivity analysis, Theory and practice in R}, SIAM, 2021.

R. Faivre, B. Iooss, S. Mahevas, D. Makowski, H. Monod, editors, 2013, \emph{Analyse de sensibilite et exploration de modeles. Applications aux modeles environnementaux}, Editions Quae.

B. Iooss and A. Saltelli, 2017, \emph{Introduction: Sensitivity analysis.} In: \emph{Springer Handbook on Uncertainty Quantification}, R. Ghanem, D. Higdon and H. Owhadi (Eds), Springer.

B. Iooss and P. Lemaitre, 2015, \emph{A review on global sensitivity analysis methods}. In \emph{Uncertainty management in Simulation-Optimization of Complex Systems: Algorithms and Applications}, C. Meloni and G. Dellino (eds), Springer.

A. Saltelli, K. Chan and E. M. Scott eds, 2000, \emph{Sensitivity Analysis}, Wiley.

A. Saltelli et al., 2008, \emph{Global Sensitivity Analysis: The Primer}, Wiley

}

\keyword{package}

