
## File: morris.R
## Description: The Morris screening method
## Author: Gilles Pujol


morris <- function(model = NULL, factors, levels, R, jump = NULL,
                   min = 0, max = 1, scale = TRUE, optim = NULL, ...){
  
  ## ARGUMENTS

  ## factor number and names
  
  if (class(factors) == "character"){
    names <- factors
    p <- length(names)
  }else{
    p <- factors
    names <- paste("X", 1:p, sep="")
  }

  ## levels

  if (length(levels) == 1){
    levels <- rep(levels, p)
  }
  
  ## boundaries
  
  if (length(min) == 1){
    min <- rep(min, p)
  }
  if (length(max) == 1){
    max <- rep(max, p)
  }
  
  ## delta

  if (is.null(jump)){
    if (1 %in% (levels %% 2)){
      stop("the levels must be even if jumps are not given...")
    }else{
      jump <- levels / 2
    }
  }else{
    if (length(jump) == 1){
      jump <- rep(jump, p)
    }
  }

  delta <- jump / (levels - 1) # corresponding to the set [0,1]^p
    
  ## DESIGN OF EXPERIMENTS

  ## B : orientation matrix
  
  B <- matrix(0, nrow = p + 1, ncol = p)
  B[lower.tri(B)] <- 1

  ## G : grid (restricted to [0, 1-delta_i])
  
  G <- rep(list(NULL), p)
  for (i in 1:p){
    G[[i]] <- (0 : (levels[i] - 1 - jump[i])) / (levels[i] - 1)
  }
  
  ## generation of the R trajectories

  x.base <- numeric(p)
  x <- NULL
  for ( i in 1 : R ){
    ## random elements : directions, and base point

    d <- sample(c(-1, 1), size = p, replace = TRUE)
    for ( i in 1 : p )
      x.base[i] <- sample(G[[i]], 1)

    ## trajectory

    B2 <- ((2 * B - 1) %*% diag(d) + 1) / 2
    x.new <- matrix(x.base, nrow = p + 1, ncol = p, byrow = TRUE) +
      0.5 * ((2 * B - 1) %*% diag(d) + 1) %*% diag(delta)
    x <- rbind(x, x.new)  
  }

  ## duplicated trajectories are removed

  x2 <- array(t(x), dim = c(p, p + 1, R))
  x2 <- unique(x2, MARGIN = 3)
  x <- matrix(x2, ncol = p, byrow = TRUE)
  R2 <- nrow(x) / (p + 1)
  if (R2 < R)
    warning(paste("keeping", R2, "trajectories out of", R))
  R <- R2

  ## optimisation of the design (selection of trajectories that optimize a
  ## good-filling-space criteria)

  if (! is.null(optim)){
    r <- optim
    x <- optim.doe(x, r)
    R <- r
  }

  ## the DOE must be transformed to be between real limits (rather than [0,1]^p)

  for ( i in 1 : ncol(x) )
    x[, i] <- x[, i] * (max[i] - min[i]) + min[i]
  colnames(x) <- names
  
  ## OBJECT OF CLASS "morris"
  
  sa <- list(model = model, levels = levels, R = R, delta = delta,
             min = min, max = max, scale = scale,
             x = x, y = NULL, ee = NULL,
             call = match.call())
  class(sa) <- "morris"
  
  ## COMPUTATION OF THE RESPONSE AND SENSITIVITY ANALYSIS
  
  if (!is.null(sa$model)){
      response(sa, ...)
      tell(sa)
  }
  
  ## RETURN OF THE OBJECT OF CLASS "morris"
 
  return(sa)
}


optim.doe <- function(x, r){
## optimisation of the design of experiment by extracting
## the r designs (r << R) that cover the space well

  p <- ncol(x)
  R <- nrow(x) / (p + 1)

  if (R < r)
    stop(paste("Can't select", r, "trajectories out of", R, sep = " "))
  
  ## indices of the trajectory #i

  traj.indices <- function(i) (1 + (i - 1) * (p + 1)) : (i * (p + 1))

  ## distances between the trajectories

  d <- matrix(0, R, R)
  for (i in 1 : (R - 1)){
    for (j in (i+1) : R){
      ## vectorized double loop...
      pi <- rep(traj.indices(i), each = p + 1)
      pj <- rep(traj.indices(j), times = p + 1)
      d[i, j] <- sum(sqrt(rowSums((x[pi, ] - x[pj, ])^2)))
    }
  }
  d <- d + t(d)

  ## selection of r trajectories

  opt <- 1
  for (i in 2 : r){
    tmp <- matrix(d[opt, -opt], nrow = length(opt))
    ## maximin criteria :
    opt.new <- which.max(apply(tmp, 2, min)) # index in (1:R)[-opt]
    opt.new <- (1:R)[-opt][opt.new] # index in 1:R
    opt <- append(opt, opt.new)
  }

  ## extraction of the matrices corresponding to the optimal trajectories

  x2 <- matrix(ncol = p, nrow = r * (p + 1))
  for (i in 1 : r){
    x2[traj.indices(i),] <- x[traj.indices(opt[i]),]
  }

  return(x2)
}


tell.morris <- function(sa, y = NULL){
  id <- deparse(substitute(sa))

  ## EXTERNAL MODEL

  if (! is.null(y))
    sa$y <- y
  
  ## SCALING OF INPUTS AND OUTPUTS (IF REQUIRED)

  x <- sa$x
  y <- sa$y
  if (sa$scale == TRUE){
    x <- as.data.frame(scale(x))
    y <- as.numeric(scale(y))
  }

  ## COMPUTING ELEMENTARY EFFECTS

  p <- ncol(x)
  ee <- matrix(nrow = sa$R, nc = p)
  colnames(ee) <- colnames(x)
  for (i in 1 : sa$R){
    j <- (1 : p) + (p + 1) * (i - 1)
    ee[i, ] <- (y[j + 1] - y[j]) / rowSums(x[j + 1, ] - x[j, ])
  }
  sa$ee <- ee
    
  ## SAVING OF THE INDICES
  
  assign(id, sa, parent.frame())
}


print.morris <- function(x, ...){
  cat("\nMORRIS OAT METHOD\n")
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  cat("\nModel runs:", length(x$y), "\n")
  if (! is.null(x$ee)){
    M <- data.frame(apply(x$ee, 2, function(x) mean(abs(x))), apply(x$ee, 2, sd))
    colnames(M) <- c("mu*", "sigma")
    rownames(M) <- colnames(x$ee)
    print(M)
  }
}


plot.morris <- function(x, identify = FALSE, ...){
  mu <- as.numeric(apply(x$ee, 2, function(x) mean(abs(x))))
  sigma <- as.numeric(apply(x$ee, 2, sd))  
  plot(mu, sigma, pch = 20,
       xlab = expression(mu^"*"), ylab = expression(sigma), ...)
  if (identify == FALSE){
    text(mu, sigma, labels = colnames(x$ee), pos = 4)
 } else{
    identify(mu, sigma, labels = colnames(x$ee))
  }
}
