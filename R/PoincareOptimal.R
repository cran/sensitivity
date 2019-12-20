# Poincare optimal constant computation for Derivative-based Global Sensitivity Measures (DGSM)           
# Authors: Olivier Roustant and Bertrand Iooss (2016)
#
# References: 
# O. Roustant, F. Barthe and B. Iooss, 
# Poincare inequalities on intervals - application to sensitivity analysis,
# Electronic Journal of Statistics, Vol. 11, No. 2, 3081-3119, 2017
#
# O. Roustant, F. Gamboa and B. Iooss, Parseval inequalities and lower bounds 
# for variance-based sensitivity indices, Preprint, ARXIV: 1906.09883

PoincareOptimal <- function (distr = list("unif", c(0, 1)), min = NULL, max = NULL, 
                             n = 500, method = c("quadrature", "integral"), only.values = TRUE, 
                             der = FALSE,  plot = FALSE, ...) {
  method <- match.arg(method, c("quadrature", "integral"))
  if ((mode(distr) != "list") & (mode(distr) != "function")) 
    stop("distr must be a list or a function")
  a <- min
  b <- max
  if (mode(distr) == "function") {
    dfct <- distr
  } else {
     if (distr[[1]] == "unif") {
      min <- distr[[2]]
      max <- distr[[3]]
      if (min >= max) stop("For the uniform distribution, 
                           the 2nd argument (min) must be smaller than the 3rd (max)")
      scaling <- max - min
      a <- 0
      b <- 1
      dfct <- dunif
    }
     if (distr[[1]] == "norm") {  # truncated normal on [min, max]
      mean <- distr[[2]]
      sd <- distr[[3]]
      if (sd <= 0) stop("For the (truncated) normal distribution, 
                        the 3rd argument (sd) must be positive")      
      a <- (min - mean)/sd
      b <- (max - mean)/sd
      scaling <- sd
      dfct <- function(x) dnorm(x) / (pnorm(b) - pnorm(a))
    }
     if (distr[[1]] == "gumbel") {   # truncated Gumbel on [min, max]
      loc <- distr[[2]]
      scale <- distr[[3]]
      if (scale <= 0) stop("For the (truncated) Gumbel distribution,
                            the 3rd argument (scale) must be positive")
      a <- (min - loc)/scale
      b <- (max - loc)/scale
      scaling <- scale
      pGumbel <- function(x) exp(-exp(-x))
      dfct <- function(x) exp(-(x + exp(-x))) / (pGumbel(b) - pGumbel(a))
    }
     if (distr[[1]] == "triangle") {
      min <- distr[[2]]
      max <- distr[[3]]
      mode <- distr[[4]]
      if (!((min < mode) & (mode < max)))
        stop("For the triangle distribution, the 2nd argument is the min, 
               the 3rd is the max, and the 4th is the mode (between min and max)")
      a <- -1
      b <- 1
      scaling <- (max - min)/2
      dfct <- function(x) (1 - abs(x)) * (abs(x) < 1)
    }
     if (distr[[1]] == "exp") {
      dfct <- function(x) {
        const <- pexp(b, distr[[2]]) - pexp(a, distr[[2]])
        dexp(x, distr[[2]]) / const
      }
    }
     if (distr[[1]] == "beta") {
      dfct <- function(x) {
        const <- pbeta(b, distr[[2]], distr[[3]]) - pbeta(a, distr[[2]], distr[[3]]) 
        dbeta(x, distr[[2]], distr[[3]]) / const
      }
    }
     if (distr[[1]] == "gamma") {
      dfct <- function(x) {
        const <- pgamma(b, distr[[2]], distr[[3]]) - pgamma(a, distr[[2]], distr[[3]])
        dgamma(x, distr[[2]], distr[[3]]) / const
      }
    }
     if (distr[[1]] == "weibull") {
      dfct <- function(x) {
        const <- pweibull(b, distr[[2]], distr[[3]]) - pweibull(a, distr[[2]], distr[[3]])
        dweibull(x, distr[[2]], distr[[3]])
      }
    }
     if (distr[[1]] == "lognorm") {
        # const <- plnorm(b, distr[[2]], distr[[3]]) - plnorm(a, distr[[2]], distr[[3]])
        dfct <- function(x) dlnorm(x, distr[[2]], distr[[3]])
     }
   }
  out <- eigenPoincare(dfct = dfct, a = a, b = b, n = n, method = method, 
                       only.values = only.values, der = der)
  if ((mode(distr) == "list") && (is.element(distr[[1]], c("unif", "triangle", "norm", "gumbel")))){
     out$opt <- out$opt * scaling^2
     out$values <- out$values / scaling^2
     out$knots <- seq(min, max, length = n)
     if (der) out$der <- out$der / scaling
     }
  if (plot & (!only.values)){
      funSol <- out$vectors[, 2]
      if (diff(funSol[1:2]) < 0) funSol <- - funSol
      plot(out$knots, funSol, type = "l", col = "blue", 
           main = paste("Solution for the Neumann problem", 
                        "\nOptimal Poincare constant :", 
                        round(out$opt, digits = 4)), 
           xlab = "t", ylab = "f(t)")
      sizeArrow <- 1/15
      arrows(x0 = min, x1 = min + (max - min) * sizeArrow, y0 = funSol[1], 
             length = sizeArrow)
      arrows(x0 = max, x1 = max - (max - min) * sizeArrow, y0 = funSol[n], 
             length = sizeArrow)
    }
  return(out)
}

# ------------------------------------------------------------

psi0 <- function(u){
  (1-abs(u))*((u>=-1) & (u<=1))
}

eigenPoincare <- function (dfct, a = 0, b = 1, n = 500, 
                           method = c("quadrature", "integral"), 
                           only.values = FALSE, der = FALSE) {
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
  
  method <- match.arg(method, c("quadrature", "integral"))
  rhoV <- dfct    # pdf of the distribution
  l <- b - a
  h <- l/n
  lower <- c(0, rep(-1, n - 1))
  upper <- c(rep(1, n - 1), 0)
  knots <- seq(a, b, length = n)
  
  if (length(method) > 1) 
    method <- "quadrature"
  ## computation of the 'rigidity' matrix
  Kh <- matrix(0, n, n)
  for (i in 1:n) {
    fUpperDiag <- function(u) {
      z <- a + (i + u) * h
      -1/h * rhoV(z)
    }
    fDiag <- function(u) {
      z <- a + (i + u) * h
      1/h * rhoV(z)
    }
    if (i < n) 
      Kh[i, i + 1] <- Kh[i + 1, i] <- integrate(fUpperDiag, 
                                                0, 1)$value
    Kh[i, i] <- integrate(fDiag, lower[i], upper[i])$value
  }
  if (method == "quadrature") {  ## quadrature with trapezoids
    ## computation of the 'mass' matrix
    W <- rhoV(seq(a, b, by = h))
    dh <- h * (W[1:n] + W[2:(n + 1)])/2
    ## resolution of the spectral problem
    dinv <- 1/sqrt(dh)
    Kh <- Kh + diag(dh, nrow = length(dh))
    A <- tcrossprod(dinv) * Kh
    spec <- eigen(A, symmetric = TRUE, only.values = only.values)
    if (!only.values) vec <- spec$vectors * dinv
  } else {
    ## computation of the 'mass' matrix
    Dh <- matrix(0, n, n)
    for (i in 1:n) {
      fUpperDiag <- function(u) {
        z <- a + (i + u) * h
        h * psi0(u) * psi0(u - 1) * rhoV(z)
      }
      fDiag <- function(u) {
        z <- a + (i + u) * h
        h * (psi0(u))^2 * rhoV(z)
      }
      if (i < n) 
        Dh[i, i + 1] <- Dh[i + 1, i] <- integrate(fUpperDiag, 
                                                  0, 1)$value
      Dh[i, i] <- integrate(fDiag, lower[i], upper[i])$value
    }
    ## resolution of the spectral problem
    L <- t(chol(Dh))
    Kh <- Kh + Dh
    Linv_Kh <- forwardsolve(L, Kh)
    A <- t(forwardsolve(L, t(Linv_Kh)))
    ## The lines above are equivalent to (but faster):
    ## Linv <- solve(L); A <- Linv%*%Kh%*%t(Linv)
    spec <- eigen(A, symmetric = TRUE, only.values = only.values)
    if (!only.values) vec <- backsolve(t(L), spec$vectors)
  }

  spec$values <- spec$values - 1
  
  # reverse order so that the first eigenvalue is 0, 
  # and the second one is the spectral gap
  spec$values <- spec$values[n:1]
  resList <- list(opt = 1 / spec$values[2], values = spec$values, knots = knots)
  if (!only.values) {
    vec <- vec[, n:1]
    resList$vectors <- vec
    if (der){
      der <- rep(0, n)
      der <- rbind(der, apply(vec, 2, diff)) / h
      dimnames(der) <- NULL
      resList$der <- der
    } 
  }
  
  return(resList)
}
