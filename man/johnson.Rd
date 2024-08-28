\name{johnson}
\alias{johnson}
\alias{print.johnson}
\alias{plot.johnson}
\alias{ggplot.johnson}

\title{Johnson indices}

\description{
  \code{johnson} computes the Johnson indices for correlated input relative importance by 
  \eqn{R^2}{R-squared} decomposition for linear and logistic regression models. These 
  indices allocates  a share of \eqn{R^2}{R-squared} to each input based on the relative 
  weight allocation (RWA) system, in the case of dependent or correlated inputs.
}

\usage{
johnson(X, y, rank = FALSE, logistic = FALSE, nboot = 0, conf = 0.95)
\method{print}{johnson}(x, \dots)
\method{plot}{johnson}(x, ylim = c(0,1), \dots)
\method{ggplot}{johnson}(data,  mapping = aes(), ylim = c(0, 1), \dots, environment
                 = parent.frame())
}

\arguments{
  \item{X}{a data frame (or object coercible by \code{as.data.frame})
    containing the design of experiments (model input variables).}
  \item{y}{a vector containing the responses corresponding to the design
    of experiments (model output variables).}
  \item{rank}{logical. If \code{TRUE}, the analysis is done on the
    ranks.}
  \item{logistic}{logical. If \code{TRUE}, the analysis is done via a
    logistic regression (binomial GLM).}
  \item{nboot}{the number of bootstrap replicates.}
  \item{conf}{the confidence level of the bootstrap confidence intervals.}
  \item{x}{the object returned by \code{johnson}.}
  \item{data}{the object returned by \code{johnson}.}
  \item{ylim}{the y-coordinate limits of the plot.}
  \item{mapping}{Default list of aesthetic mappings to use for plot. If not specified, 
    must be supplied in each layer added to the plot.}
  \item{environment}{[Deprecated] Used prior to tidy evaluation.}
  \item{\dots}{arguments to be passed to methods, such as graphical
    parameters (see \code{par}).}
}

\value{
  \code{johnson} returns a list of class \code{"johnson"}, containing the following
  components:

  \item{call}{the matched call.}
  \item{johnson}{a data frame containing the estimations of the johnson
    indices, bias and confidence intervals.}
}

\details{
  Logistic regression model (\code{logistic = TRUE}) and rank-based indices
  (\code{rank = TRUE}) are incompatible.
}

\references{

  L. Clouvel, B. Iooss, V. Chabridon, M. Il Idrissi and F. Robin, 2024,
  \emph{An overview of variance-based importance measures in the linear regression context: 
  comparative analyses and numerical tests}, Preprint.
  \url{https://hal.science/hal-04102053}
  
B. Iooss, V. Chabridon and V. Thouvenot, \emph{Variance-based importance 
measures for machine learning model interpretability}, Congres lambda-mu23,
Saclay, France, 10-13 octobre 2022
\url{https://hal.science/hal-03741384}

J.W. Johnson, 2000, \emph{A heuristic method for estimating the relative 
weight of predictor variables in multiple regression}, Multivariate 
Behavioral Research, 35:1-19.

J.W. Johnson and J.M. LeBreton, 2004, \emph{History and use of relative 
importance indices in organizational research}, Organizational 
Research Methods, 7:238-257.

}

\author{
Bertrand Iooss and Laura Clouvel
}

\examples{

##################################
# Same example than the one in src()

# a 100-sample with X1 ~ U(0.5, 1.5)
#                   X2 ~ U(1.5, 4.5)
#                   X3 ~ U(4.5, 13.5)

library(boot)
n <- 100
X <- data.frame(X1 = runif(n, 0.5, 1.5),
                X2 = runif(n, 1.5, 4.5),
                X3 = runif(n, 4.5, 13.5))

# linear model : Y = X1 + X2 + X3

y <- with(X, X1 + X2 + X3)

# sensitivity analysis

x <- johnson(X, y, nboot = 100)
print(x)
plot(x)

library(ggplot2)
ggplot(x)


#################################
# Same examples than the ones in lmg()

library(boot)
library(mvtnorm)

set.seed(1234)
n <- 1000
beta<-c(1,-1,0.5)
sigma<-matrix(c(1,0,0,
                0,1,-0.8,
                0,-0.8,1),
              nrow=3,
              ncol=3)

##########
# Gaussian correlated inputs

X <-rmvnorm(n, rep(0,3), sigma)
colnames(X)<-c("X1","X2", "X3")

#########
# Linear Model

y <- X\%*\%beta + rnorm(n,0,2)

# Without Bootstrap confidence intervals
x<-johnson(X, y)
print(x)
plot(x)

# With Boostrap confidence intervals
x<-johnson(X, y, nboot=100, conf=0.95)
print(x)
plot(x)

# Rank-based analysis
x<-johnson(X, y, rank=TRUE, nboot=100, conf=0.95)
print(x)
plot(x)

#######
# Logistic Regression
y<-as.numeric(X\%*\%beta + rnorm(n)>0)
x<-johnson(X,y, logistic = TRUE)
plot(x)
print(x)

#################################
# Test on a modified Linkletter fct with: 
# - multivariate normal inputs (all multicollinear)
# - in dimension 50 (there are 42 dummy inputs)
# - large-size sample (1e4)

library(mvtnorm)

n <- 1e4
d <- 50
sigma <- matrix(0.5,ncol=d,nrow=d)
diag(sigma) <- 1
X <- rmvnorm(n, rep(0,d), sigma)

y <- linkletter.fun(X)
joh <- johnson(X,y)
sum(joh$johnson) # gives the R2
plot(joh)
}

\seealso{
\code{\link{src}}, \code{\link{lmg}}, \code{\link{pmvd}}, \code{\link{johnsonshap}}
}

\keyword{regression}
