\name{sensitivity}
\alias{sensitivity-package}
\alias{sensitivity}
\docType{package}

\title{Package 'sensitivity': Sensitivity Analysis}

\description{
  The \pkg{sensitivity} package implements sensitivity analysis methods: linear
  and monotonic sensitivity analysis (SRC, PCC, SRRC, PRCC), the
  screening method of Morris, and non-linear global sensitivity analysis
  (the Sobol indices, the FAST method).

  The functions of this package generate the design of experiments
  (depending on the method of analysis) and compute the sensitivity
  indices based on the model inputs and outputs. All sensitivity indices
  can be estimated with the bootstrap technique which allows to estimate
  the bias, and basic bootstrap confidence intervals. Text and graphical
  outputs display the results of the analysis.
}

\details{
  The approach applied when performing a sensitivity analysis (SA) is as
  follows:
  \describe{
    \item{step 1}{The model is defined: it is a function that associates
      an unique real value (the response), with input parameters (the
      factors).}
    \item{step 2}{A sensitivity analysis method is chosen. Parameters of
      this method must be in accordance with objectives and technical
      constraints (like computational time).}
    \item{step 3}{A design of experiments (DOE) corresponding to the SA
      method is generated.}
    \item{step 4}{The model is evaluated on the DOE values.}
    \item{step 5}{The sensitivity indices are computed, based on input
      and output values.}
    \item{step 6}{Post-treatments\ldots}
  }

  The \pkg{sensitivity} package allows to follow this methodology:

  (step 1) The model can be internal or external to R. If internal, it
  can be a function that takes an unique \code{matrix} or \code{data.frame}
  parameter and returns a \code{numeric} vector. It can also be a
  predictor, i.e. an object wich can be called with the \code{predict}
  method. If the model is external it does not have to be
  interfaced with R: the program will stop after generating the DOE, and
  calculations will start again when one gives the corresponding
  responses.
  
  The four next steps depend upon the type of the model:

  For internal models:\cr
  (step 2-5) \code{sa <- method(model, parameters\dots)}
  
  For external models:\cr
  (step 2-3) \code{sa <- method(model = NULL, parameters\dots)}\cr
  (step 4) external to R, and the result is loaded in the y variable\cr
  (step 5) \code{compute(sa, y)}


  \code{method} should be the name of a SA function, such as
  \code{\link{linsa}}, \code{\link{morris}}, \code{\link{sobol}},
  \code{\link{sobol.sal02}} or \code{\link{fast}}. These function create
  the object \code{sa} of class \code{"linsa"}, \code{"morris"},
  \code{"sobol"}, \code{"sobol.sal02"} or \code{"fast"}. For further
  information on these function, see the corresponding documentation.  

  Finally, for displaying the results of the analysis:\cr
  (step 6) \code{print(sa); plot(sa)}
}

\references{
Saltelli, A., Chan, K. and Scott, E. M., 2000, \emph{Sensitivity
  analysis}. Wiley, 33--45.
}

\seealso{
  \code{\link{linsa}}
  \code{\link{morris}}
  \code{\link{sobol}}  
  \code{\link{sobol.sal02}}
  \code{\link{fast}}
  \code{\link{compute}}
  \code{\link{testmodels}}
}

\keyword{misc}

