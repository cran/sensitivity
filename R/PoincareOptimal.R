# Poincare optimal constant computation for Derivative-based Global Sensitivity Measures (DGSM)           
# Authors: Olivier Roustant and Bertrand Iooss (2016)
#
# Reference: 
# O. Roustant, F. Barthe and B. Iooss, 
# Poincare inequalities on intervals - application to sensitivity analysis,
# Submitted, 2016 - https://hal.archives-ouvertes.fr/hal-01388758

PoincareOptimal <- function(distr=list("unif",c(0,1)), min = NULL, max = NULL, n = 500, 
                            method = c("quadrature", "integral"), only.values = TRUE,
                            plot = FALSE, ...){
  
  if( (mode(distr) != "list") & (mode(distr) != "function")) stop("distr must be a list or a function")
  
  if (mode(distr) == "function"){
    dfct <- distr
    a = min ; b = max
  } else{
    
    if (distr[[1]] == "unif"){
      dfct <- function(x) rep(1, length(x))
      a = 0 ; b = 1 ; scaling = (distr[[3]]-distr[[2]])^2
    }
    
    if (distr[[1]] == "norm"){
      dfct <- dnorm
      mean = distr[[2]] ; sd = distr[[3]]
      a = (min-mean)/sd ; b = (max-mean)/sd ; scaling = sd^2
    }
    
    if (distr[[1]] == "gumbel"){
      dfct = function(x) exp(- ( x + exp(-x) ) )
      loc = distr[[2]] ; scale = distr[[3]]
      a = (min - loc)/scale ; b = (max - loc)/scale ; scaling = scale^2
    }
    if (distr[[1]] == "triangle"){
      dab = function(x){
        cond <- (x < distr[[4]])
        cond * (x-distr[[2]])/(distr[[4]]-distr[[2]]) + (!cond) * (distr[[3]]-x)/(distr[[3]]-distr[[4]])
      }
      dfct = function(x) dab(distr[[2]] + x * (distr[[3]]-distr[[2]]))
      a = 0 ; b = 1 ; scaling = (distr[[3]]-distr[[2]])^2
    }
    
    if (distr[[1]] == "exp"){
      dfct <- function(x) dexp(x,distr[[2]])
      a = min ; b = max
    }
    
    if (distr[[1]] == "beta"){
      dfct <- function(x) dbeta(x,distr[[2]],distr[[3]])
      a = min ; b = max
    }
    
    if (distr[[1]] == "gamma"){
      dfct = function(x) dgamma(x,distr[[2]],distr[[3]])
      a = min ; b = max
    }
    
    if (distr[[1]] == "weibull"){
      dfct = function(x) dweibull(x,distr[[2]],distr[[3]])
      a = min ; b = max
    }
    
    if (distr[[1]] == "lognorm"){
      dfct = function(x) dlnorm(x,distr[[2]],distr[[3]])
      a = min ; b = max
    }
  }
  
  out <- eigenPoincare(dfct=dfct, a=a, b=b, n=n, method=method, only.values=only.values, plot=plot)
  
  if ((mode(distr) == "list") && (is.element(distr[[1]],c("unif","norm","triangle","gumbel")))){
    out$opt = out$opt * scaling
    out$values = out$values * scaling 
  }

  return(out) # list(opt = 1/spec$values[n-1], values = spec$values, vectors = vec)
                                                
}

# ------------------------------------------------------------

psi0 <- function(u){
  (1-abs(u))*((u>=-1) & (u<=1))
}

eigenPoincare <- function(dfct, a = 0, b = 1, n = 500, method = c("quadrature", "integral"), 
                          only.values = FALSE, plot = TRUE){
## This function solves numerically the spectral problem corresponding to
## the Poincare inequality, with Neumann conditions.
## The differential equation is  f'' - V'f'  = - lambda f  
## with f'(a) = f'(b) = 0.
## ---
## V : potential, actually such that the pdf dfct = exp(-V), 
#     up to a mutliplicative constant (which has no impact on computations)
# method : "integral" is longer but does not do any approximation
#          "quadrature" uses the trapez quadrature (close and quicker)
# only.values : to compute only the eigen values
# plot : if TRUE plots a minimizer of the Rayleigh ratio

rhoV <- dfct   # pdf of the distribution
l <- b - a 
h <- l / n

lower <- c(0, rep(-1, n-1))
upper <- c(rep(1, n-1), 0)

if (length(method) > 1) method <- "quadrature"

## computation of the 'rigidity' matrix
Kh <- matrix(0, n, n)
for (i in 1:n){
  fUpperDiag <- function(u){
    z <- a + (i + u) * h
    - 1/h * rhoV(z)
  }
  fDiag <- function(u){
    z <- a + (i + u) * h
    1/h * rhoV(z)
  }
  if (i < n) Kh[i, i+1] <- Kh[i+1, i] <- integrate(fUpperDiag, 0, 1)$value
  Kh[i, i] <- integrate(fDiag, lower[i], upper[i])$value
}

if (method == "quadrature"){  ## quadrature with trapezoids
  ## computation of the 'mass' matrix
  W <- rhoV(seq(a, b, by = h))
  dh <- h * (W[1:n] + W[2:(n+1)]) / 2
  ## resolution of the spectral problem
  dinv <- 1 / sqrt(dh)
  Kh <- Kh + diag(dh, nrow = length(dh))
  A <- tcrossprod(dinv) * Kh
  spec <- eigen(A, symmetric=TRUE, only.values = only.values)
  vec <- spec$vectors * dinv
} else {
  ## computation of the 'mass' matrix
  Dh <- matrix(0, n, n)
  for (i in 1:n){
    fUpperDiag <- function(u){
      z <- a + (i + u) * h
      h * psi0(u) * psi0(u-1) * rhoV(z)
    }
    fDiag <- function(u){
      z <- a + (i + u) * h
      h * (psi0(u))^2  * rhoV(z)
    }
    if (i < n) Dh[i, i+1] <- Dh[i+1, i] <- integrate(fUpperDiag, 0, 1)$value
    Dh[i, i] <- integrate(fDiag, lower[i], upper[i])$value
  }
  ## resolution of the spectral problem
  L <- t(chol(Dh))
  Kh <- Kh + Dh
  Linv_Kh <- forwardsolve(L, Kh)
  A <- t(forwardsolve(L, t(Linv_Kh)))
  ## The lines above are equivalent to (but faster):
  #Linv <- solve(L); A <- Linv%*%Kh%*%t(Linv)
  spec <- eigen(A, symmetric=TRUE, only.values = only.values)
  if (!only.values){
    vec <- backsolve(t(L), spec$vectors)
  }
  else{
    vec <- NULL
  }
}

spec$values <- spec$values - 1

normalize <- function(x) x / sqrt(crossprod(x))
if (!only.values) vec <- apply(vec, 2, normalize)

## plot of the eigen vector corresponding to minimal eigen value
if (plot & (!only.values)) {
  funSol <- vec[, n-1]
  if (diff(funSol[1:2]) < 0) funSol <- - funSol   ## there exists a non-negative solution
  plot(seq(a, b, length = n), funSol, 
       type = "l", col = "blue",
       main = paste("Solution for the Neumann problem",
                    "\nOptimal Poincare constant :", 
                    round(1/spec$values[n-1], digits = 4)),
       xlab = "t", ylab = "f(t)")
  sizeArrow <- 1/15
  arrows(x0 = a, x1 = a + (b-a)*sizeArrow, y0 = funSol[1], 
         length = sizeArrow)
  arrows(x0 = b, x1 = b - (b-a)*sizeArrow, y0 = funSol[n], 
         length = sizeArrow)
}

return(list(opt = 1/spec$values[n-1], values = spec$values, vectors = vec))
}
