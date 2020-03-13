\name{squaredIntEstim}
\alias{squaredIntEstim}

\title{
Squared integral estimate
}
\description{
This function provides two estimators of a squared expectation. 
The first one, naive, is the square of the sample mean. It is positively biased.
The second one is a U-statistics, and unbiased. The two are equivalent for large sample sizes.
}
\usage{
squaredIntEstim(x, method = "unbiased")
}
\arguments{
  \item{x}{A vector of observations supposed to be drawn independently from a square integrable random variable}
  \item{method}{If "unbiased", computes the U-statistics, otherwise the square of the sample mean is computed}
}

\details{
Let X1, ..., Xn be i.i.d. random variables. The aim is to estimate t = E(Xi)^2.
The naive estimator is the square of the sample mean: T1 = [(X1 + ... + Xn)/n]^2.
It is positively biased, and the bias is equal to s^2/n, where s^2 = var(X1).
The U-statistics estimator is the average of Xi * Xj over all unordered pairs (i,j).
Equivalently, it is equal to T1 minus the (unbiased) sample variance divided by n.
}

\value{
A real number, corresponding to the estimated value of the squared integral.
}

\references{
O. Roustant, F. Gamboa and B. Iooss, \emph{Parseval inequalities and lower bounds for 
variance-based sensitivity indices}, Electronic Journal of Statistics, 14:386-412, 2020

Van der Vaart, A. W. Asymptotic statistics. Vol. 3. Cambridge university press, 2000.
}

\author{
O. Roustant
}

\examples{
n <- 100  # sample size
nsim <- 100 # number of simulations
mu <- 0

T <- Tunb <- rep(NA, nsim)
theta <- mu^2  # E(X)^2, with X following N(mu, 1)

for (i in 1:nsim){
  x <- rnorm(n, mean = mu, sd = 1)
  T[i] <- squaredIntEstim(x, method = "biased")
  Tunb[i] <- squaredIntEstim(x, method = "unbiased")
}

par(mfrow = c(1, 1))
boxplot(cbind(T, Tunb))
abline(h = theta, col = "red")
abline(h = c(mean(T), mean(Tunb)), col = c("blue", "cyan"), lty = "dotted")
# look at the difference between median and mean
}