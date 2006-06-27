\name{linsa}
\alias{linsa}
\alias{compute.linsa}
\alias{print.linsa}
\alias{plot.linsa}

\title{Linear Sensitivity Analysis}

\description{
  \code{linsa} computes the standardized regression coefficients (SRC)
  and the partial correlation coefficients (PCC). Analysis can be done
  on the ranks; then the indices are the standardized rank regression
  coefficients (SRRC) and the partial rank correlation coefficients
  (PRCC).  
}

\usage{
linsa(model = NULL, x, pcc = TRUE, rank = FALSE,
      nboot = 0, conf = 0.95, \dots)
\method{compute}{linsa}(sa, y = NULL)
}

\arguments{
  \item{model}{the model.}
  \item{x}{the input sample.}
  \item{pcc}{logical. If \code{TRUE}, the P(R)CCs are computed.}
  \item{rank}{logical. If \code{TRUE}, the analysis is done on the
    ranks.}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level for bootstrap confidence intervals.}
  \item{sa}{the sensitivity analysis object.}
  \item{y}{the response.}
  \item{\dots}{any other arguments for \code{model} which are passed
    unchanged each time it is called.}
}

\details{
  \code{model} is a function or a predictor (a class with a
  \code{predict} method) computing the response \code{y} based on the
  sample given by \code{x}. If no model is specified, the indices will be
  computed when one gives the response.
}

\value{
  \code{linsa} returns an object of class \code{"linsa"}.
  An object of class \code{"linsa"} is a list containing the following
  components:

  \item{model}{the model.}
  \item{x}{the input sample.}
  \item{rank}{logical. If \code{TRUE}, the analysis was done on the
    ranks.}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level for bootstrap confidence intervals.}
  \item{y}{the response.}
  \item{src}{the estimations of the SRC indices (or SRRC if
    rank analysis is requested).}
  \item{pcc}{if requested, the estimations of the PCC
    indices (or PRCC if rank analysis is requested).}
  \item{call}{the matched call.}
}

\references{
Saltelli, A., Chan, K. and Scott, E. M., 2000, \emph{Sensitivity
  analysis}, Wiley.
}

\seealso{
  \code{\link{sensitivity}}
  \code{\link{compute}}
}

\examples{
# linear model : Y = X1 + X2 + X3

model1 <- function(x) x[, 1] + x[, 2] + x[, 3]

# a 500-sample with X1 ~ U(0.5, 1.5)
#                   X2 ~ U(1.5, 4.5)
#                   X3 ~ U(4.5, 13.5)

n <- 500
x <- data.frame(X1 = runif(n, 0.5, 1.5),
                X2 = runif(n, 1.5, 4.5),
                X3 = runif(n, 4.5, 13.5))

# sensitivity analysis

sa <- linsa(model = model1, x = x, nboot = 100)
print(sa)
par(mfrow = c(1,2))
plot(sa, ask = FALSE)
}

\keyword{misc}
