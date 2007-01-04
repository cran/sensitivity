
## File: fast_saltelli99.R
## Description: FAST method, as improved by Saltelli (so-called "extended-FAST")
## Author: Gilles Pujol


fast.saltelli99 <- function(model = NULL, factors, n, M = 4, omega = NULL,
                            q = NULL, q.arg = NULL, ...){
  ## ARGUMENTS

  ## factor number and names

  if (class(factors) == "character"){
    names <- factors
    p <- length(names)
  }
  else{
    p <- factors
    names <- paste("X", 1:p, sep="")
  }

  ## quantiles

  if (is.null(q)){
    q <- "qunif"
  }
  if (is.null(q.arg)){
    q.arg <- list()
  }
  if (length(q) == 1){
    q <- rep(q, p)
    q.arg <- rep(list(q.arg), p)
  }

  ## DESIGN OF EXPERIMENTS

  ## set of frequencies

  if (is.null(omega)){
    omega <- numeric(p)
    omega[1] <- floor((n-1) / (2 * M))
    m <- floor(omega[1] / (2 * M))
    if (m >= p - 1){
      omega[-1] <- floor(seq(from = 1, to = m, length.out = p - 1))
    }else{
      omega[-1] <- 0 : (p - 2) %% m + 1
    }
  }

  ## discretization of the s-space
  
  s <- 2 * pi / n * (0 : (n - 1))

  ## transformation to get points in the x-space
  
  x <- matrix(nrow = n * p, ncol = p)
  omega2 <- numeric(p)
  for (i in 1 : p){
    omega2[i] <- omega[1]
    omega2[-i] <- omega[-1]
    l <- (1 : n) + (i - 1) * n
    for (j in 1 : p){
      g <- 0.5 + 1 / pi * asin(sin(omega2[j] * s))
      x[l, j] <- do.call(q[j], c(list(p = g), q.arg[[j]]))
    }
  }
  colnames(x) <- names

  ## OBJECT OF CLASS "fast.saltelli99"

  sa <- list(method = "saltelli99", model = model, M = M,
             s = s, omega = omega, x = data.frame(x), y = NULL, D = NULL,
             D1 = NULL, Dt = NULL, call = match.call())
  class(sa) <- c("fast.saltelli99")
  
  ## COMPUTATION OF THE RESPONSE AND SENSITIVITY ANALYSIS
  
  if (!is.null(sa$model)){
      response(sa, ...)
      tell(sa)
  }
  
  ## RETURN OF THE OBJECT OF CLASS "fast"
  
  return(sa)
}


estim.fast.saltelli99 <- function(y, omega, M){
  n <- nrow(y)
  p <- ncol(y) - 1 # (remark : s is in column number p+1)
  D <- numeric(p)
  D1 <- numeric(p)
  Dt <- numeric(p)
  for (j in 1 : p){
    f <- fft(y[, j], inverse = FALSE)
    s <- (Mod(f[2 : (n / 2)]) / n)^2
    D[j] <- 2 * sum(s)
    D1[j] <- 2 * sum(s[(1 : M) * omega[1]])
    Dt[j] <- 2 * sum(s[1 : (omega[1] / 2)]) # [attention, source de debordement potentielle...]
  }
  return(c(D, D1, Dt))
}


tell.fast.saltelli99 <- function(sa, y = NULL){
  id <- deparse(substitute(sa))

  ## EXTERNAL MODEL

  if (! is.null(y))
    sa$y <- y
  
  ## ESTIMATION OF THE INDICES

  p <- ncol(sa$x)
  n <- length(sa$s)
  
  data <- cbind(matrix(sa$y, nr = n, nc = p, byrow = FALSE), sa$s)
  indices <- estim.fast.saltelli99(data, sa$omega, sa$M)
  sa$D <- indices[1 : p]
  sa$D1 <- indices[(p + 1) : (2 * p)]
  sa$Dt <- indices[(2 * p + 1) : (3 * p)]

  ## SAVING OF THE INDICES
  
  assign(id, sa, parent.frame())
}


print.fast.saltelli99 <- function(x, ...){
  cat("\nFOURIER AMPLITUDE SENSITIVITY TEST\n")
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  cat("\nModel runs:", length(x$y), "\n")
  M <- cbind(x$D1 / x$D, 1 - x$Dt / x$D)
  colnames(M) <- c("first order", "total order")
  rownames(M) <- colnames(x$x)
  cat("\nEstimations of the indices:\n")
  print(M)
}


plot.fast.saltelli99 <- function(x, ...){
  M <- rbind(x$D1 / x$D, 1 - x$Dt / x$D - x$D1 / x$D)
  colnames(M) <- colnames(x$x)
  barplot(M, ylim=c(0,1), col=c("white","grey"))
  legend("topright", c("single effect", "interactions"), fill=c("white","grey"))
}
