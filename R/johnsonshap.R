# Johnson-Shapley indices
#
# Bertrand Iooss 2023

estim.johnsonshap <- function(model, data, N, i = 1:nrow(data)){
    dd <- data[i, ]
    d <- ncol(dd)

    # Computation of the weight matrix
    cor_matrix <- cor(dd, use = "pairwise.complete.obs")
    corXX <- cor_matrix # Correlations between inputs
    corXX.eigen <- eigen(corXX) # Eigenvalues and eigenvectors
    W <- corXX.eigen$vec %*% sqrt(diag(corXX.eigen$val)) %*% t(corXX.eigen$vec)
    
    # Computation of indices (alpha) between output and orthogonal variables
    # With the orthogonalized inputs Z
    PDQ <- svd(dd)
    P <- PDQ$u
    Q <- PDQ$v
    Z <- P %*% t(Q)
    
    # on calcule les indices de Sobol de Z (en les transformant en X)
    sobrepZ <- sensitivity::sobolrep(model=NULL, factors=d, N=N)
    for (i in 1:d) sobrepZ$X[,i] <- quantile(ecdf(Z[,i]), sobrepZ$X[,i])
    sobrepX <- sobrepZ$X %*% W
    colnames(sobrepX) <- colnames(sobrepZ$X)
    yrepZ <- model(sobrepX)
    tell(sobrepZ, yrepZ)
    # effets de Shapley (= 1er ordre + somme 2eme ordre / 2)
    alpha <- sobrepZ$S$original # 1er ordre
    a <- gtools::combinations(d,2,1:d) # on cherche les 2nd ordre
    for (i in 1:d){
      nr <- which(a == i)%%nrow(a)
      nr[nr==0] <- nrow(a)
      alpha[i] <- alpha[i] + sum(sobrepZ$S2$original[nr]) / 2
    }
    #alpha <- alpha / sum(alpha) # normalisation pas tres catholique
    rsquare <- sum(alpha)
    alpha <- matrix(alpha)
    
    
#    W1 <- Q %*% diag(PDQ$d) %*% t(Q)
#    W[,1:ncol(W1)] <-  W1[,1:ncol(W1)] / sqrt(rowSums(W1^2))[1:ncol(W1)]
    
    W^2 %*% alpha ^ 2
    
}


johnsonshap <- function(model = NULL, X, N, rank = FALSE, nboot = 0, conf = 0.95) {

  if (rank) {
    for (i in 1:ncol(X)) {
      data[,i] <- rank(X[,i])
    }
  }
  
  if (nboot == 0) {
    johnsonshap <- data.frame(original = estim.johnsonshap(model, X, N))
    rownames(johnsonshap) <- colnames(X)
  } else {
    boot.johnsonshap <- boot(X, estim.johnsonshap, R = nboot, model = model, N = N)
    johnsonshap <- bootstats(boot.johnsonshap, conf, "basic")
    rownames(johnsonshap) <- colnames(X)
  }

  out <- list(model = model, X = X, N = N, rank = rank, nboot = nboot, conf = conf,
              call = match.call())
  class(out) <- "johnsonshap"
  out$johnsonshap <- johnsonshap
    
  return(out)
}


print.johnsonshap <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  cat("\nJohnson-Shapley indices:\n")
  print(x$johnsonshap)
}


plot.johnsonshap <- function(x, ylim = c(0,1), ...) {  
  nodeplot(x$johnsonshap, ylim = ylim, main = "Johnson-Shapley indices")
}

ggplot.johnsonshap <- function(data, mapping = aes(), ylim = c(0,1), ..., environment = parent.frame()) {  
  x <- data
  nodeggplot(listx = list(x$johnsonshap), xname = "Johnson-Shapley indices", ylim = ylim, title = "Johnson-Shapley indices")
}
