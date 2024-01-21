\name{EPtest}
\alias{EPtest}

\title{Non-parametric variable significance test based on the empirical process}

\description{
  \code{EPtest} builds the non-parametric variable significance test from Klein and Rochet (2022) for the null hypothesis \eqn{H_0: S^u = S}{H0: Su = S} where \eqn{S^u}{Su} is the Sobol index for the inputs \eqn{X_i, i \in u}{Xi,i in u} ans \eqn{S}{S} is the Sobol index for all the inputs in \eqn{X}{X}. 
}

\usage{
EPtest(X, y, u = NULL, doe = NULL, Kdoe = 10, tau = 0.1)
}

\arguments{
  \item{X}{a matrix or data.frame that contains the numerical inputs as columns.}
  \item{y}{a vector of output.}
	\item{u}{the vector of indices of the columns of X for which we want to test the significance.}
	\item{doe}{the design of experiment on which the empirical process is to be evaluated. It should be independent from X.}
	\item{Kdoe}{if doe is null and Kdoe is specified, the design of experiment is taken as Kdoe points drawn uniformly independently on intervals delimited by the range of each input.}
	\item{tau}{a regularization parameter to approximate the limit chi2 distribution of the test statistics under H0.}
}

\value{
  \code{EPtest} returns a list containing:
  \item{statistics}{The test statistics that follows a chi-squared distribution under the null hypothesis.}
  \item{ddl}{The number of degrees of freedom used in the limit chi-square distribution for the test.}
  \item{p-value}{The test p-value.}
}

\references{
T. Klein and P. Rochet, \emph{Test comparison for Sobol Indices over nested sets of variables}, SIAM/ASA Journal on Uncertainty Quantification 10.4 (2022): 1586-1600.}

\author{Paul Rochet}

\seealso{\link{sobol}}

\examples{

# Model: Ishigami
  
n = 100
X = matrix(runif(3*n, -pi, pi), ncol = 3)
  
y = ishigami.fun(X)
	
# Test the significance of X1, H0: S1 = 0
EPtest(X[, 1], y, u = NULL)

# Test if X1 is sufficient to explain Y, H0: S1 = S123
EPtest(X, y, u = 1)
  
# Test if X3 is significant in presence of X2, H0: S2 = S23
EPtest(X[, 2:3], y, u = 1)
  
}
