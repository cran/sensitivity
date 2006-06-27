fast <- function(model = NULL, factors, n, M = 4,
                 G = "uniform", min = 0, max = 1,
                 nboot = 0, conf = 0.95, ...)
{
  # ARGUMENTS

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

  omega <- numeric(p)
  omega[1] <- floor((n-1) / (2 * M))
  m <- floor(omega[1] / (2 * M))
  if (m >= p - 1)
    omega[-1] <- floor(seq(from = 1, to = m, length.out = p - 1))
  else
    omega[-1] <- 0 : (p - 2) %% m + 1

  # frequencies for each index
  
  Omega <- matrix(nrow = p, ncol = p)
  for (i in 1 : p){
    Omega[i, i] <- omega[1]
    Omega[i, -i] <- omega[-1]
  }

  # uniform sampling of the s-space
  
  s <- runif(n, min = -pi, max = pi)

  # transformation to get points in the x-space
  
  x <- matrix(nrow = n * p, ncol = p)
  for (i in 1 : p){
    for (j in 1 : p){
      l <- (1 : n) + (i - 1) * n
      x[l, j] <- G(i, sin(Omega[i, j] * s))
    }
  }
  colnames(x) <- names

  # OBJECT OF CLASS "fast"

  sa <- list(model = model, M = M, nboot = nboot, conf = conf,
             s = s, omega = Omega, x = data.frame(x), y = NULL, S = NULL,
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


estim.extendedfast <- function(data, i, omega, M)
{
  d <- data[i, ]
  n <- nrow(d)
  p <- ncol(d) - 1 # (remark : s is in column number p+1)
  S <- numeric(p)
  St <- numeric(p)
  for (j in 1 : p){
    theta = d[, p+1] %o% as.numeric((1 : M) %o% omega[j, ])
    y = matrix(rep(d[, j], M * p), n, M * p)
    A = matrix(1 / n * colSums(y * cos(theta)), M, p)
    B = matrix(1 / n * colSums(y * sin(theta)), M, p)
    D = 2 * sum(A^2 + B^2)
    S[j] = 2 * sum(A[, j]^2 + B[, j]^2) / D
    St[j] = 1 - 2 * sum(A[, -j]^2 + B[, -j]^2) / D
  }
  return(c(S, St))
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
  indices <- estim(estim.extendedfast, data, sa$nboot, conf = sa$conf,
                   omega = sa$omega, M = sa$M)

  sa$S <- data.frame(indices[1:p,])
  colnames(sa$S) <- colnames(indices)
  rownames(sa$S) <- paste("S", 1:p, sep = "")

  sa$St <- data.frame(indices[(p+1):(2*p),])
  colnames(sa$St) <- colnames(indices)
  rownames(sa$St) <- paste("St", 1:p, sep = "")

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
  p <- ncol(sa$x)
  nodeplot(x$S, xlim = c(1,p+1), ylim = c(0,1), labels = colnames(x$x),
           pch = pch[1])
  nodeplot(x$St, xlim = c(1,p+1), ylim = c(0,1), labels = FALSE,
           pch = pch[2], at = (1:p)+.3, add = TRUE)
  legend(x = "topright", legend = c("first order", "total order"), pch = pch)
}
