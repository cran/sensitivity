\name{PoincareChaosSqCoef}
\alias{PoincareChaosSqCoef}

\title{Squared coefficients computation in generalized chaos}

\description{
  This program computes the squared coefficient of the function decomposition in the tensor basis formed by eigenfunctions of Poincare differential operators.
  After division by the variance of the model output, it provides lower bounds of first-order and total Sobol' indices.
}

\usage{
PoincareChaosSqCoef(PoincareEigen, multiIndex, design, output, outputGrad = NULL, 
                    inputIndex = 1, der = FALSE, method = "unbiased")
}

\arguments{
\item{PoincareEigen}{output list from PoincareOptimal() function}
\item{multiIndex}{vector of indices (l1, ..., ld). A coordinate equal to 0 corresponds to the constant basis function 1.}
\item{design}{design of experiments (matrix of size n x d) with d the number of inputs and n the number of observations}
\item{output}{vector of length n (y1, ..., yn) of output values at \code{design} points}
\item{outputGrad}{matrix n x d whose columns contain the output partial derivatives at \code{design} points}
\item{inputIndex}{index of the input variable (between 1 and d)}
\item{der}{logical (default=FALSE): should we use the formula with derivatives to compute the squared coefficient?}
\item{method}{"biased" or "unbiased" formula when estimating the squared integral. See \code{\link{squaredIntEstim}}}
}

\details{
Similarly to polynomial chaos, where tensors of polynomials are used, we consider here tensor 
basis formed by eigenfunctions of Poincare differential operators. This basis is also orthonormal, 
and Parseval formula lead to lower bound for (unnormalized) Sobol, total Sobol indices, and any variance-based index.
Denoting by \eqn{(e_{1, l1}... e_{d, ld})} one tensor basis, the corresponding coefficient is equal to

\eqn{c_{l1, ..., ld} = <f, e_{1, l1}... e_{d, ld}>}.

For a given input variable (say \eqn{x1} to simplify notations), it can be rewritten with derivatives as:

\eqn{c_{l1, ..., ld} = <df/dx1, de_{1, l1}/dx1 e_{2, l2}...e_{d, ld}>.}

The function returns an estimate of \eqn{c_{l1, ..., ld}^2}, corresponding to one of these two forms (derivative-free, or derivative-based).
}

\value{
An estimate of the squared coefficient.
}

\references{
O Roustant, F. Gamboa, B Iooss. Parseval inequalities and lower bounds 
for variance-based sensitivity indices. 2019. hal-02140127
}

\author{
Olivier Roustant and Bertrand Iooss
}

\seealso{
  \code{\link{PoincareOptimal}}
}

\examples{

# A simple example

g <- function(x, a){
  res <- x[, 1] + a*x[, 1]*x[, 2]
  attr(res, "grad") <- cbind(1 + a * x[, 2], a * x[, 1])
  return(res)
}

n <- 1e3
set.seed(0)
X <- matrix(runif(2*n, min = -1/2, max = 1/2), nrow = n, ncol = 2)
a <- 3
fX <- g(X, a = a)

out_1 <- out_2 <- PoincareOptimal(distr = list("unif", -1/2, 1/2), 
                                  only.values = FALSE, der = TRUE, 
                                  method = "quad")
out <- list(out_1, out_2)

# Lower bounds for X1
c2_10 <- PoincareChaosSqCoef(PoincareEigen = out, multiIndex = c(1, 0), 
                             design = X, output = fX, outputGrad = attr(fX, "grad"), 
                             inputIndex = 1, der = FALSE)
c2_11 <- PoincareChaosSqCoef(PoincareEigen = out, multiIndex = c(1, 1), 
                             design = X, output = fX, outputGrad = attr(fX, "grad"), 
                             inputIndex = 1, der = FALSE)
c2_10_der <- PoincareChaosSqCoef(PoincareEigen = out, multiIndex = c(1, 0), 
                                 design = X, output = fX, outputGrad = attr(fX, "grad"), 
                                 inputIndex = 1, der = TRUE)
c2_11_der <- PoincareChaosSqCoef(PoincareEigen = out, multiIndex = c(1, 1), 
                                 design = X, output = fX, outputGrad = attr(fX, "grad"), 
                                 inputIndex = 1, der = TRUE)

LB1 <- c(8/pi^4, c2_10, c2_10_der)
LB1tot <- LB1 + c(64/pi^8 * a^2, c2_11, c2_11_der)
LB <- cbind(LB1, LB1tot)
rownames(LB) <- c("True lower bound value", 
                  "Estimated, no derivatives", "Estimated, with derivatives")
colnames(LB) <- c("D1", "D1tot")
cat("True values of D1 and D1tot:", c(1/12, 1/12 + a^2 / 144),"\n")
cat("Sample size: ", n, "\n")
cat("Lower bounds computed with the first Poincare eigenvalue:\n")
print(LB)
cat("\nN.B. Increase the sample size to see the convergence to true lower bound values.\n")

############################################################
# Flood model example (see Roustant et al., 2017, 2019)

\donttest{

library(evd) # Gumbel law
library(triangle) # Triangular law

# Flood model
Fcrues_full2=function(X,ans=0){ 
  # ans=1 gives Overflow output; ans=2 gives Cost output; ans=0 gives both
  mat=matrix(X,ncol=8);
  if (ans==0){ reponse=matrix(NA,nrow(mat),2);}
  else{ reponse=rep(NA,nrow(mat));}
  for (i in 1:nrow(mat)) {
    H = (mat[i,1] / (mat[i,2]*mat[i,8]*sqrt((mat[i,4] - mat[i,3])/mat[i,7])))^(0.6) ;
    S = mat[i,3] + H - mat[i,5] - mat[i,6] ;
    if (S > 0){ Cp = 1 ;}
    else{ Cp = 0.2 + 0.8 * (1 - exp(-1000 / S^4));}
    if (mat[i,5]>8){ Cp = Cp + mat[i,5]/20 ;}
    else{ Cp = Cp + 8/20 ;}
    if (ans==0){
      reponse[i,1] = S ;
      reponse[i,2] = Cp ;
    }
    if (ans==1){ reponse[i] = S ;}
    if (ans==2){ reponse[i] = Cp ;}
    
  }
  return(RES=reponse)
}

# Flood model derivatives (by finite-differences)
dFcrues_full2 <- function(X, i, ans, eps){
  der = X
  X1 = X
  X1[,i] = X[,i]+eps
  der = (Fcrues_full2(X1,ans) - Fcrues_full2(X,ans))/(eps)
  return(der)
}

# Function for flood model inputs sampling
EchantFcrues_full2<-function(taille){
  X = matrix(NA,taille,8)
  X[,1] = rgumbel.trunc(taille,loc=1013.0,scale=558.0,min=500,max=3000)
  X[,2] = rnorm.trunc(taille,mean=30.0,sd=8,min=15.)
  X[,3] = rtriangle(taille,a=49,b=51,c=50)
  X[,4] = rtriangle(taille,a=54,b=56,c=55)
  X[,5] = runif(taille,min=7,max=9)
  X[,6] = rtriangle(taille,a=55,b=56,c=55.5)
  X[,7] = rtriangle(taille,a=4990,b=5010,c=5000)
  X[,8] = rtriangle(taille,a=295,b=305,c=300)
  return(X)
}

d <- 8
n <- 1e3
eps <- 1e-7 # finite-differences for derivatives
x <- EchantFcrues_full2(n)
yy <- Fcrues_full2(x, ans=2)
y <- scale(yy, center = TRUE, scale = FALSE)[,1]
dy <- NULL
for (i in 1:d) dy <- cbind(dy, dFcrues_full2(x, i, ans=2, eps))

method <- "quad"
out_1 <- PoincareOptimal(distr = list("gumbel", 1013, 558), min=500,max=3000, 
                         only.values = FALSE, der = TRUE, method = method)
out_2 <- PoincareOptimal(distr = list("norm", 30, 8), min=15, max=200, 
                         only.values = FALSE, der = TRUE, method = method)
out_3 <- PoincareOptimal(distr = list("triangle", 49, 51, 50), 
                         only.values = FALSE, der = TRUE, method = method)
out_4 <- PoincareOptimal(distr = list("triangle", 54, 56, 55), 
                         only.values = FALSE, der = TRUE, method = method)
out_5 <- PoincareOptimal(distr = list("unif", 7, 9), 
                         only.values = FALSE, der = TRUE, method = method)
out_6 <- PoincareOptimal(distr = list("triangle", 55, 56, 55.5), 
                         only.values = FALSE, der = TRUE, method = method)
out_7 <- PoincareOptimal(distr = list("triangle", 4990, 5010, 5000), 
                         only.values = FALSE, der = TRUE, method = method)
out_8 <- PoincareOptimal(distr = list("triangle", 295, 305, 300), 
                         only.values = FALSE, der = TRUE, method = method)
out_ <- list(out_1,out_2,out_3,out_4,out_5,out_6,out_7,out_8)

c2 <- c2der <- c2tot <- c2totder <- rep(0,d)

for (i in 1:d){
  m <- diag(1,d,d) ; m[,i] <- 1
  
  for (j in 1:d){
    cc <- PoincareChaosSqCoef(PoincareEigen = out_, multiIndex = m[j,], 
            design = x, output = y, outputGrad = NULL, 
            inputIndex = i, der = FALSE)
    c2tot[i] <- c2tot[i] + cc
    if (j == i) c2[i] <- cc
    
    cc <- PoincareChaosSqCoef(PoincareEigen = out_, multiIndex = m[j,], 
            design = x, output = y, outputGrad = dy, 
            inputIndex = i, der = TRUE)
    c2totder[i] <- c2totder[i] + cc
    if (j == i) c2der[i] <- cc
  }
}

print("Lower bounds of first-order Sobol' indices without derivatives:")
print(c2/var(y))
print("Lower bounds of first-order Sobol' indices with derivatives:")
print(c2der/var(y))

print("Lower bounds of total Sobol' indices without derivatives:")
print(c2tot/var(y))
print("Lower bounds of total Sobol' indices with derivatives:")
print(c2totder/var(y))

}

}
