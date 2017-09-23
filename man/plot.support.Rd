\name{plot.support}
\alias{plot.support}
\alias{plot}
\alias{scatterplot.support}
\alias{scatterplot}

\title{
Support index functions: Measuring the effect of input variables over their support
}

\description{
Methods to plot the normalized support index functions (Fruth et al., 2016). 
}

\usage{
\method{plot}{support}(x, i = 1:ncol(x$X),
        xprob = FALSE, p = NULL, p.arg = NULL,
        ylim = NULL, col = 1:3, lty = 1:3, lwd = c(2,2,1), cex = 1, ...)
\method{scatterplot}{support}(x, i = 1:ncol(x$X), 
               xprob = FALSE, p = NULL, p.arg = NULL, 
               cex = 1, cex.lab = 1, ...)

}

\arguments{
  \item{x}{an object of class support.}
  \item{i}{an optional vector of integers indicating the subset of input variables \code{X_i} for plotting. Default is the entire set of input variables.}
  \item{xprob}{an optional boolean indicating whether the inputs should be plotted in probability scale.}
  \item{p}{,}
  \item{p.arg}{list of probability names and parameters for the input distribution.}
  \item{ylim}{,}
  \item{col}{,}
  \item{lty}{,}
  \item{lwd}{,}
  \item{cex}{,}
  \item{cex.lab}{usual graphical parameters.}
  \item{\dots}{additional graphical parameters to be passed to \code{scatterplot} method (\code{ggMarginal} function).}
}

\details{
If \code{xprob = TRUE}, the input variable \code{X_i} is plotted in probability scale according to the informations provided in the arguments \code{p, p.arg}: The x-axis is thus \code{F(x)}, where \code{F} is the cdf of \code{X_i}. If these ones are not provided, the empirical distribution is used for rescaling: The x-axis is thus \code{Fn(x)}, where \code{Fn} is the empirical cdf of \code{X_i}.

Legend details:

zeta*T : normalized total support index function

zeta* : normalized 1st-order support index function

nu* : normalized DGSM

Notice that the sum of (normalized) DGSM (nu*) over all input variables is equal to 1. 
Furthermore, the expectation of the total support index function (zeta*T) is equal to the (normalized) DGSM (nu*).
}

\author{O. Roustant}

\seealso{ 
Estimation of support index functions: \code{\link{support}} 
}
