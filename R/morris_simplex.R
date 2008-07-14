# Morris' simplex sub-routines. (See main file morris.R)
#
# Gilles Pujol 2007-2008


simplex.reg <- function(p) {
# generates the matrix of a regular simplex, of edge length = 1,
# centered on the origin
  S <- matrix(0, nr = p + 1, nc = p)
  S[2,1] <- 1
  for (i in 3 : (p + 1)) {
    for (j in 1 : (i - 2)) {
      S[i,j] <- mean(S[1 : (i - 1), j])
    }
    S[i, i - 1] <- sqrt(1 - sum(S[i, 1 : (i - 2)]^2))
  }
  scale(S, scale = FALSE)
}


simplex.rect <- function(p) {
# generates the matrix of a rectangle ("orthonormal") simplex
  S <- matrix(0, nr = p + 1, nc = p)
  for (i in 1:p) {
    S[i+1,i] <- 1
  }
  S
}


plane.rot <- function(p, i, j, theta) {
# matrix of the plane (i,j)-rotation of angle theta in dimension p
  R <- diag(nr = p)
  R[c(i,j), c(i,j)] <- matrix(c(cos(theta), sin(theta), - sin(theta), cos(theta)), nr = 2)
  return(R) 
}

random.simplexes <- function(p, r, min = rep(0, p), max = rep(1, p), h = 0.25) {
# generates r random simplexes

  X <- matrix(nrow = r * (p + 1), ncol = p)
  
  # initial random rotation matrix
  R <- diag(nr = p, nc = p)
  ind <- combn(p, 2)
  theta <- runif(choose(p, 2), min = 0, max = 2 * pi)
  for (i in 1 : choose(p, 2)) {
    R <- plane.rot(p, ind[1,i], ind[2,i], theta[i]) %*% R
  }
  
  # reference simplex
  S.ref <- R %*% (h * t(simplex.reg(p)))
  
  for (i in 1 : r) {
    # one more random plane rotation of the reference simplex
    ind <- sample(p, 2)
    theta <- runif(1, min = 0, max = 2 * pi)
    S.ref <- plane.rot(p, ind[1], ind[2], theta) %*% S.ref
    
    # translation to put the simplex into the domain
    tau.min <- min - apply(S.ref, 1, min)
    tau.max <- max - apply(S.ref, 1, max)
    X[ind.rep(i, p),] <- t(S.ref + runif(n = p, min = tau.min, max = tau.max))
  }

  return(X)
}


ee.simplex <- function(X, y) {
# compute the elementary effects for a simplex design
  p <- ncol(X)
  r <- nrow(X) / (p + 1)
  ee <- matrix(nr = r, nc = p)
  colnames(ee) <- colnames(X)
  for (i in 1 : r) {
    j <- ind.rep(i, p)
    ee[i, ] <- solve(cbind(as.matrix(X[j,]), rep(1, p + 1)), y[j])[1 : p]
  }
  return(ee)
}
