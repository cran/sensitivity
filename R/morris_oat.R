# Morris' OAT sub-routines. (See main file morris.R)
#
# Gilles Pujol 2006


random.oat <- function(p, r, binf = rep(0, p), bsup = rep(0, p), nl, design.step) {
    # orientation matrix B
    B <- matrix(-1, nrow = p + 1, ncol = p)
    B[lower.tri(B)] <- 1
    # grid step
    delta <- design.step / (nl - 1)
    X <- matrix(nrow = r * (p + 1), ncol = p)
    for (j in 1 : r) {
        # directions matrix D
        D <- diag(sample(c(-1, 1), size = p, replace = TRUE))
        # permutation matrix P
        perm <- sample(p)
        P <- matrix(0, nrow = p, ncol = p)
        for (i in 1 : p) {
            P[i, perm[i]] <- 1
        }
        # starting point
        x.base <- matrix(nrow = p + 1, ncol = p)
        for (i in 1 : p) {
            x.base[,i] <- ((sample(nl[i] - design.step[i], size = 1) - 1) / (nl[i] - 1))
        }
        X[ind.rep(j,p),] <- 0.5 * (B %*% P %*% D + 1) %*% diag(delta) + x.base
    }
    for (i in 1 : p) {
        X[,i] <- X[,i] * (bsup[i] - binf[i]) + binf[i]
    }
    return(X)
}

ee.oat <- function(X, y) {
# compute the elementary effects for a OAT design
  p <- ncol(X)
  r <- nrow(X) / (p + 1)
  ee <- matrix(nrow = r, ncol = p)
  colnames(ee) <- colnames(X)
  for (i in 1 : r) {
    j <- ind.rep(i, p)
    j1 <- j[1 : p]
    j2 <- j[2 : (p + 1)]
    #ee[i,] <- (y[j2] - y[j1]) / rowSums(X[j2,] - X[j1,])
    ee[i,] <- solve(X[j2,] - X[j1,], y[j2] - y[j1])
  }
  return(ee)
}
