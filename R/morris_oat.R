# Morris' OAT sub-routines. (See main file morris.R)
#
# Gilles Pujol 2006


random.oat <- function(p, nl, jump) {
  delta <- jump / (nl - 1) 
  B <- matrix(0, nrow = p + 1, ncol = p)                 # orientation matrix
  B[lower.tri(B)] <- 1
  x.base <- matrix(nrow = p + 1, ncol = p)               # base point
  D <- diag(sample(c(-1, 1), size = p, replace = TRUE))  # directions 
  for (i in 1 : p) {                                    
    x.base[,i] <- ((sample(nl[i] - jump[i], size = 1) - 1) / (nl[i] - 1))
  }
  0.5 * ((2 * B - 1) %*% D + 1) %*% diag(delta) + x.base
}


ee.oat <- function(X, y) {
# compute the elementary effects for a OAT design
  p <- ncol(X)
  r <- nrow(X) / (p + 1)
  ee <- matrix(nr = r, nc = p)
  colnames(ee) <- colnames(X)
  for (i in 1 : r) {
    j <- ind.rep(i, p)
    j1 <- j[1 : p]
    j2 <- j[2 : (p + 1)]
    ee[i,] <- (y[j2] - y[j1]) / rowSums(X[j2,] - X[j1,])
  }
  return(ee)
}
