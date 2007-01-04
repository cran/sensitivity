\name{srcpcc}
\alias{srcpcc}
\alias{tell.srcpcc}
\alias{print.srcpcc}
\alias{plot.srcpcc}

\title{Linear Sensitivity Analysis}

\description{
  \code{srcpcc} computes the standardized regression coefficients (SRC)
  and the partial correlation coefficients (PCC). Analysis can be done
  on the ranks; then the indices are the standardized rank regression
  coefficients (SRRC) and the partial rank correlation coefficients
  (PRCC).  
}

\usage{
srcpcc(model = NULL, x, pcc = TRUE, rank = FALSE,
      nboot = 0, conf = 0.95, \dots)
}

\arguments{
  \item{model}{the model}
  \item{x}{the input sample}
  \item{pcc}{logical. If \code{TRUE}, the P(R)CCs are computed}
  \item{rank}{logical. If \code{TRUE}, the analysis is done on the
    ranks}
  \item{nboot}{the number of bootstrap replicates}
  \item{conf}{the confidence level for bootstrap confidence intervals}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called}
}

\details{
  \code{model} is a function or a predictor (a class with a
  \code{predict} method) computing the response \code{y} based on the
  sample given by \code{x}. If no model is specified, the indices will be
  computed when one gives the response.
}

\value{
  \code{srcpcc} returns an object of class \code{"srcpcc"}.
  An object of class \code{"srcpcc"} is a list containing the following
  components:

  \item{y}{the response}
  \item{src}{the estimations of the SRC indices (or SRRC if
    rank analysis is requested)}
  \item{pcc}{if requested, the estimations of the PCC
    indices (or PRCC if rank analysis is requested)}
}

\section{Computational cost}{
  The number of model evaluations is \eqn{n}{n} where \eqn{n}{n} is the
  size of the sample \code{x}.
}

\references{
Saltelli, A., Chan, K. and Scott, E. M., 2000, \emph{Sensitivity
  analysis}, Wiley.
}

\examples{
# linear model : Y = X1 + X2 + X3

model1 <- function(x) x[, 1] + x[, 2] + x[, 3]

# a 100-sample with X1 ~ U(0.5, 1.5)
#                   X2 ~ U(1.5, 4.5)
#                   X3 ~ U(4.5, 13.5)

n <- 100
x <- data.frame(X1 = runif(n, 0.5, 1.5),
                X2 = runif(n, 1.5, 4.5),
                X3 = runif(n, 4.5, 13.5))

# sensitivity analysis

sa <- srcpcc(model = model1, x = x, nboot = 100)
print(sa)
par(mfrow = c(1,2))
plot(sa, ask = FALSE)
}

\keyword{misc}
