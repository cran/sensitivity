# Morris's screening method (Morris 1992, Campolongo 2007)
# Provides also simplex-based screening designs (Pujol 2007)
#
# Gilles Pujol 2006-2008
#
# Sub-files:
# * morris_oat.R
# * morris_simplex.R
# * morris_sfd.R


ind.rep <- function(i, p) {
# indices of the ith trajectory in the DoE
  (1 : (p + 1)) + (i - 1) * (p + 1)
}


morris <- function(model = NULL, factors, r, design, binf = 0, bsup = 1,
                   scale = TRUE, ...) {
  
  # argument checking: factor number and names
  if (is.character(factors)) {
    X.labels <- factors
    p <- length(X.labels)
  } else if (is.numeric(factors)) {
    p <- factors
    X.labels <- paste("X", 1 : p, sep="")
  } else {
    stop("invalid argument \'factors\', waiting for a scalar (number) or a character string vector (names)")
  }

  # argument checking: number of repetitions
  if (length(r) == 1) {
    r.max <- r
  } else {
    r.max <- r[2]
    r <- r[1]
  }

  # argument checking: design parameters
  if (! "type" %in% names(design)) {
    design$type <- "oat"
    warning("argument \'design$type\' not found, set at \'oat\'")
  }
  if (design$type == "oat") {
    # one-at-a-time design
    if (! "levels" %in% names(design)) {
      stop("argument \'design$levels\' not found")
    }
    nl <- design$levels
    if (length(nl) == 1) nl <- rep(nl, p)
    if ("grid.jump" %in% names(design)) {
      jump <- design$grid.jump
      if (length(jump) == 1) jump <- rep(jump, p)
    } else {
      jump <- rep(1, p)
      warning("argument \'design$grid.jump\' not found, set at 1")
    }
  } else if (design$type == "simplex") {
    # simplex-based design
    h <- design$scale.factor
  } else {
    stop("invalid argument design$type, waiting for \"oat\" or \"simplex\"")
  }
    
  # argument checking: domain boundaries
  if (length(binf) == 1) binf <- rep(binf, p)
  if (length(bsup) == 1) bsup <- rep(bsup, p)
  
  # generation of the initial design
  if (design$type == "oat") {
    X <- matrix(nrow = r.max * (p + 1), ncol = p)
    for (j in 1 : r.max) {
      X[ind.rep(j,p),] <- random.oat(p, nl, jump)
    }
    for (i in 1 : p) {
      X[,i] <- X[,i] * (bsup[i] - binf[i]) + binf[i]
    }
  } else if (design$type == "simplex") {
    X <- random.simplexes(p, r.max, binf, bsup, h)
  }
  
  # duplicated repetitions are removed
  X.unique <- array(t(X), dim = c(p, p + 1, r.max))
  X.unique <- unique(X.unique, MARGIN = 3)
  X <- matrix(X.unique, ncol = p, byrow = TRUE)
  colnames(X) <- X.labels
  r.unique <- nrow(X) / (p + 1)
  if (r.unique < r.max) {
    warning(paste("keeping", r.unique, "repetitions out of", r.max))
  }
  r.max <- r.unique

  # optimization of the design
  if (r < r.max) {
    ind <- morris.maximin(X, r)
    X <- X[sapply(ind, function(i) ind.rep(i, p)),]
  }

  # object of class "morris"
  x <- list(model = model, factors = factors, r = r, design = design,
            binf = binf, bsup = bsup, scale = scale, X = X, call =
            match.call())
  class(x) <- "morris"

  # computing the response if the model is given
  if (!is.null(x$model)) {
    response(x, ...)
    tell(x)
  }
   
  return(x)
}


tell.morris <- function(x, y = NULL, ...) {
  id <- deparse(substitute(x))

  if (! is.null(y)) {
    x$y <- y
  } else if (is.null(x$y)) {
    stop("y not found")
  }

  X <- x$X
  y <- x$y

  if (x$scale) {
    X <- scale(X)
    y <- as.numeric(scale(y))
  }

  if (x$design$type == "oat") {
    x$ee <- ee.oat(X, y)
  } else if (x$design$type == "simplex") {
    x$ee <- ee.simplex(X, y)
  }

  assign(id, x, parent.frame())
}


print.morris <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (! is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    mu <- apply(x$ee, 2, mean)
    mu.star <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma <- apply(x$ee, 2, sd)
    
    out <- data.frame(mu, mu.star, sigma)
    rownames(out) <- colnames(x$ee)
    print(out)
  } else {
    cat("(empty)\n")
  }
}


plot.morris <- function(x, identify = FALSE, ...) {
  if (! is.null(x$ee)) {
    mu.star <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma <- apply(x$ee, 2, sd)
    
    plot(mu.star, sigma, pch = 20, xlab = expression(mu^"*"),
         ylab = expression(sigma), ...)
  
    if (identify) {
      identify(mu.star, sigma, labels = colnames(x$ee))
    } else {
      text(mu.star, sigma, labels = colnames(x$ee), pos = 4)
    }
  }
}


plot3d.morris <- function(x, alpha = c(0.2, 0), sphere.size = 1) {
  library(rgl)
  spheres.rad <- max((apply(x$ee,2,max) - apply(x$ee,2,min)) / 100) * sphere.size
  color = "grey"
  cone.nfaces = 100
  
  mu <- apply(x$ee, 2, mean)
  mu.star <- apply(x$ee, 2, function(x) mean(abs(x)))
  sigma <- apply(x$ee, 2, sd)
    
  open3d()
  
  xmax <- max(mu.star)
  zmax <- max(max(sigma), xmax)
  
  n <- 100
  theta <- seq(from = 0, to = pi, length.out = n + 1)
  x <- rep(c(0, xmax, xmax), n)
  y <- as.numeric(rbind(rep(0, n), - xmax * cos(theta[-(n + 1)]),
                        - xmax * cos(theta[-1])))
  z <- as.numeric(rbind(rep(0, n), xmax * sin(theta[-(n + 1)]),
                        xmax * sin(theta[-1])))
  triangles3d(x, y, z, color = color, alpha = alpha[1])
  
  x <- rep(c(0, xmax, xmax, 0), 2)
  y <- c(0, -xmax, -xmax, 0, 0, xmax, xmax, 0)
  z <- c(0, 0, zmax, zmax, 0, 0, zmax, zmax)
  quads3d(x, y, z, color = color, alpha = alpha[2])
  
  plot3d(mu.star, mu, sigma, type = "s", radius = spheres.rad, add = TRUE)
  
  axes3d()
  title3d(xlab="mu*", ylab="mu", zlab="sigma")
}
