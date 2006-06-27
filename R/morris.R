morris <- function(model = NULL, factors, levels, r, k.delta = "usual",
                   min = 0, max = 1, scale = TRUE,
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

  # delta

  if (k.delta == "usual")
    if (p %% 2 == 0)
      k.delta <- p / 2

  delta <- k.delta / (p - 1)
    
  # DESIGN OF EXPERIMENTS

  matrices <- function(p, levels, delta)
  {
    # D.star : p*p diagonal matrix composed of equiprobable +1 and -1
    
    D.star <- diag(2 * floor(runif(p, 0, 2)) - 1)
   
    # P.star : p*p random permutation matrix
    
    tmp <- sample(p)
    P.star <- matrix(0, p, p)
    for ( j in 1:p )
      P.star[j, tmp[j]] <- 1
    
    # B : (p+1)*p matrix with ones in the lower triangular part and zeros in the
    # upper part
    
    B <- matrix(0, p, p)
    B[lower.tri(B, diag = TRUE)] <- 1
    B <- rbind(matrix(0, 1, p), B)
    
    # B.star : (p+1)*p matrix of the DOE
    
    M <- delta / 2 * ((2 * B - matrix(1, p + 1, p)) %*%
                      D.star + matrix(1, p + 1, p))
    B.star <- 1.1
    while ( max(B.star) > 1 ){
      x <- floor(runif(p, 0, levels - 1)) / (levels - 1)
      B.star <- (matrix(1, p + 1, 1) %*% x + M) %*% P.star
    }
    
    return(list(B = B.star, D = D.star, P = P.star))
  }

  # generation of all matrices (r repetitions)
  
  B <- matrix(nr = (p + 1) * r, nc = p)
  D <- matrix(nr = p * r, nc = p)
  P <- matrix(nr = p * r, nc = p)
  for ( i in 1:r ){
    tmp <- matrices(p, levels, delta)
    B[(1 : (p + 1)) + (i - 1) * (p + 1),] <- tmp$B
    D[(1 : p) + (i - 1) * p,] <- tmp$D
    P[(1 : p) + (i - 1) * p,] <- tmp$P
  }

  # the DOE must be transformed to be between real limits (rather than [0,1]^p)

  x <- B
  for (i in 1:ncol(x))
    x[,i] <- B[,i] * (max[i] - min[i]) + min[i]
  colnames(x) <- names
  
  # OBJECT OF CLASS "morris"
  
  sa <- list(model = model, levels = levels, r = r, delta = delta,
             min = min, max = max, scale = scale, nboot = nboot, conf = conf,
             B = B, D = D, P = P,  x = x, y = NULL, mu = NULL, sigma = NULL,
             call = match.call())
  class(sa) <- "morris"
  
  # COMPUTATION OF THE RESPONSE AND SENSITIVITY ANALYSIS
  
  if (!is.null(sa$model)){
      response(sa, ...)
      compute(sa)
  }
  
  # RETURN OF THE OBJECT OF CLASS "morris"
 
  return(sa)
}


estim.morris <- function(data, i) # data : (r*p) matrix of the
                                  # elementary effects
{
  d <- data[i, ]
  return(as.numeric(c(mean(abs(d)), sd(d))))
}


compute.morris <- function(sa, y = NULL)
{
  id <- deparse(substitute(sa))

  # EXTERNAL MODEL

  if (! is.null(y))
    sa$y <- y
  
  # SCALING OF INPUTS AND OUTPUTS (IF REQUIRED)
  
  if (sa$scale == TRUE){
    sa$x <- as.data.frame(scale(sa$x))
    sa$y <- as.numeric(scale(sa$y))
  }
  
  # delta vector (computed on the factor sample)

  p <- ncol(sa$x)
  down <- (1 : (sa$r * (p + 1)))[- ((0 : (sa$r - 1)) * (p + 1) + 1)]
  up <- (1 : (sa$r * (p + 1)))[- ((0 : (sa$r - 1)) * (p + 1) + (p + 1))]
  delta <- as.numeric(abs(t(sa$x[up, ] - sa$x[down, ])[t(sa$P) == 1]))

  # d : sample of the elementary effects
  
  d <- matrix(nr=sa$r, nc=p)
  colnames(d) <- colnames(sa$x)
  for ( i in 1:sa$r ){
    y <- sa$y[(1 : (p + 1)) + (i - 1) * (p + 1)]
    D <- sa$D[(1 : p) + (i - 1) * p,]
    P <- sa$P[(1 : p) + (i - 1) * p,]
    d[i, ] <- (((y[2 : (p + 1)] - y[1 : p]) /
                      delta[(1 : p) + (i - 1) * p]) %*% D %*% P)
  }
  d <- as.data.frame(d)
  
  # ESTIMATION OF THE INDICES
  
  estimation <- estim(estim.morris, d, sa$nboot, sa$conf)

  sa$mu <- as.data.frame(estimation[1:p,])
  colnames(sa$mu) <- colnames(estimation)
  rownames(sa$mu) <- colnames(sa$x)

  sa$sigma <- as.data.frame(estimation[(p+1):(2*p),])
  colnames(sa$sigma) <- colnames(estimation)
  rownames(sa$sigma) <- colnames(sa$x)
  
  # SAVING OF THE INDICES
  
  assign(id, sa, parent.frame())
}


print.morris <- function(x, ...)
{
  cat("\nMORRIS OAT METHOD\n")
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  cat("\nModel runs:", length(x$y), "\n")
  cat("\nmu*:\n")
  print(x$mu)
  cat("\nsigma:\n")
  print(x$sigma)
}


plot.morris <- function(x, ...)
{
  xlab <- expression(mu^"*")
  ylab <- expression(sigma)
  crossplot(x$mu, x$sigma, col = "lightgray", xlab = xlab, ylab = ylab,
            labels = colnames(x$x))
}

