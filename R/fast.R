fast <- function(model = NULL, factors, n, M = 4,
                 G = "uniform", min = 0, max = 1, omega = NULL, ...)
{
  # ARGUMENTS..

  # factor number and names

  if (class(factors) == "character"){
    names <- factors
    p <- length(names)
  }
  else{
    p <- factors
    names <- paste("X", 1:p, sep="")
  }

  # boundaries
  
  if (length(min) == 1)
    min <- rep(min, p)
  if (length(max) == 1)
    max <- rep(max, p)
  
  # space transformation for uniform factors

  if (G == "uniform"){
    G <- function(i, x)
      a[i] + (b[i] - a[i]) * (0.5 + 1 / pi * asin(x))
    environment(G) <- new.env()
    assign("a", min, environment(G))
    assign("b", max, environment(G))
  }
  
  # DESIGN OF EXPERIMENTS

  # set of frequencies

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

  # discretization of the s-space
  
  #s <- runif(n, min = -pi, max = pi)
  s <- 2 * pi / n * (0 : (n - 1))

  # transformation to get points in the x-space
  
  x <- matrix(nrow = n * p, ncol = p)
  omega2 <- numeric(p)
  for (i in 1 : p){
    omega2[i] <- omega[1]
    omega2[-i] <- omega[-1]
    l <- (1 : n) + (i - 1) * n
    for (j in 1 : p){
      x[l, j] <- G(j, sin(omega2[j] * s))
    }
  }
  colnames(x) <- names

  # OBJECT OF CLASS "fast"

  sa <- list(model = model, M = M,
             s = s, omega = omega, x = data.frame(x), y = NULL, S = NULL,
             St = NULL, call = match.call())
  class(sa) <- "fast"
  
  # COMPUTATION OF THE RESPONSE AND SENSITIVITY ANALYSIS
  
  if (!is.null(sa$model)){
      response(sa, ...)
      compute(sa)
  }
  
  # RETURN OF THE OBJECT OF CLASS "fast"
  
  return(sa)
}


## estim.fast <- function(data, i, omega, M)
## {
##   d <- data[i, ]
##   n <- nrow(d)
##   p <- ncol(d) - 1 # (remark : s is in column number p+1)
##   s <- d[, p + 1]
##   S <- numeric(p)
##   St <- numeric(p)
##   for (j in 1 : p){
##     theta <-  omega[j,] %o% (1 : M)
##     y <- d[, j]
##     A <- matrix(colSums(y * cos(s %o% as.numeric(theta))) / n, nr = p, nc = M)
##     B <- matrix(colSums(y * sin(s %o% as.numeric(theta))) / n, nr = p, nc = M)
##     D <- var(y) # OUAIS, BON....
##     S[j] <- 2 * sum(A[j,]^2 + B[j,]^2) / D
##     St[j] <- 1 - 2 * sum(A[-j,]^2 + B[-j,]^2) / D # FAUX
##   }
##   return(c(S, St))
## }


estim.fastfft <- function(y, omega, M){
  n <- nrow(y)
  p <- ncol(y) - 1 # (remark : s is in column number p+1)
  D <- numeric(p)
  D1 <- numeric(p)
  Dt <- numeric(p)
  for (j in 1 : p){
    F <- fft(y[, j], inverse = FALSE)
    T <- length(F)
    Spectre <- Mod(F[2 : (T /2)])^2 / T^2
    print(Spectre)
    D[j] <- 2 * sum(Spectre)
    D1[j] <- 2 * sum(Spectre[(1 : M) * omega[1]])
    Dt[j] <- 2 * sum(Spectre[1 : (omega[1] / 2)]) # attention, source de debordement potentielle...
  }
  print(D)
  print(D1)
  print(Dt)
  return(c(D1 / D, 1 - Dt / D))
}


compute.fast <- function(sa, y = NULL)
{
  id <- deparse(substitute(sa))

  # EXTERNAL MODEL

  if (! is.null(y))
    sa$y <- y
  
  # ESTIMATION OF THE INDICES

  p <- ncol(sa$x)
  n <- length(sa$s)
  
  data <- cbind(matrix(sa$y, nr = n, nc = p, byrow = FALSE), sa$s)
##   indices <- estim(estim.fast, data, sa$nboot, conf = sa$conf,
##                    omega = sa$omega, M = sa$M)
  indices <- data.frame(original = estim.fastfft(data, sa$omega, sa$M))

  sa$S <- data.frame(indices[1:p,])
  colnames(sa$S) <- colnames(indices)
##   rownames(sa$S) <- paste("S", 1:p, sep = "")
  rownames(sa$S) <- names(sa$x)

  sa$St <- data.frame(indices[(p+1):(2*p),])
  colnames(sa$St) <- colnames(indices)
##   rownames(sa$St) <- paste("St", 1:p, sep = "")
  rownames(sa$St) <- names(sa$x)

  # SAVING OF THE INDICES
  
  assign(id, sa, parent.frame())
}


print.fast <- function(x, ...){
  cat("\nFOURIER AMPLITUDE SENSITIVITY TEST\n")
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  cat("\nModel runs:", length(x$y), "\n")
  cat("\nFirst order indices:\n")
  print(x$S)
  cat("\nTotal order indices:\n")
  print(x$St)
}


plot.fast <- function(x, ...){
  pch <- c(21, 24)
  p <- ncol(x$x)
  nodeplot(x$S, xlim = c(1,p+1), ylim = c(0,1), labels = colnames(x$x),
           pch = pch[1])
  nodeplot(x$St, xlim = c(1,p+1), ylim = c(0,1), labels = FALSE,
           pch = pch[2], at = (1:p)+.3, add = TRUE)
  legend(x = "topright", legend = c("first order", "total order"), pch = pch)
}
