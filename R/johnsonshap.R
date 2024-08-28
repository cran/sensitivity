# Johnson-Shapley indices
#
# Bertrand Iooss 2024

johnsonshap <- function(model = NULL, X1, N, nboot = 0, conf = 0.95) {
    d <- dim(X1)[[2]]

    # Computation of indices (alpha) between output and orthogonal variables
    # With the orthogonalized inputs Z
    PDQ <- svd(X1)
    P <- PDQ$u # left singular vectors
    Q <- PDQ$v # right singular vectors
    Z <- P %*% t(Q)
    
    # Computation of the weight matrix
    D <- PDQ$d # singular values
    W <- Q %*% diag(D) %*% t(Q)
    W2sum <- sqrt(colSums(W^2))
    Wstar <- t(t(W) / W2sum)
    
    # computing Sobol' indices of Z (by transformation of Z to X)
    sobrepZ <- sensitivity::sobolrep(model=NULL, factors=d, N=N, nboot = nboot)
    for (i in 1:d) sobrepZ$X[,i] <- quantile(ecdf(Z[,i]), sobrepZ$X[,i])
    sobrepX <- sobrepZ$X %*% W
    colnames(sobrepX) <- colnames(sobrepZ$X)
    X <- sobrepX
    
    x <- list(model = model, X1 = X1, N = N, nboot = nboot, conf = conf,
              X = X, sobrepZ = sobrepZ, Wstar = Wstar, call = match.call())
    class(x) <- "johnsonshap"
    
    #calcul of the response for explicit model
    if (! is.null(x$model)){
      response(x)
      tell(x, x$y)
    }  
    return(x)
}

tell.johnsonshap <- function(x, y = NULL, ...){
  
    id <- deparse(substitute(x))
    if (! is.null(y)) {
      x$y <- y
    } 
    else if (is.null(x$y)) {
      stop("y not found")
    }
    
    sob <- x$sobrepZ
    sensitivity::tell(sob, x$y)
    
    # Shapley effects(= 1st order + sum of 2nd order / 2)
    d <- dim(x$X1)[[2]]
    a <- gtools::combinations(d,2,1:d) # look for 2nd order indices

    if (x$nboot == 0) {
      alpha <- sob$S$original # 1er ordre
      for (i in 1:d){
        nr <- which(a == i)%%nrow(a)
        nr[nr==0] <- nrow(a)
        alpha[i] <- alpha[i] + sum(sob$S2$original[nr]) / 2
      }
      alpha <- matrix(alpha)
      estim <- x$Wstar^2 %*% alpha
    
      johnsonshap <- data.frame(original = estim, row.names = colnames(x$X1))
    } else {
      estim <- matrix(0, nrow = d, ncol = 5)
      for (k in c(1,4,5)){
        alpha <- sob$S[,k] # 1st order
        for (i in 1:d){
          nr <- which(a == i)%%nrow(a)
          nr[nr==0] <- nrow(a)
          alpha[i] <- alpha[i] + sum(sob$S2$original[nr]) / 2
        }
        alpha <- matrix(alpha)
        estim[,k] <- x$Wstar^2 %*% alpha
      }
      johnsonshap <- data.frame(estim, row.names = colnames(x$X1))
      colnames(johnsonshap) <-  c("original", "bias", "std. error", "min. c.i.", "max. c.i.")
    }
    
    x$sobrepZ <- sob
    x$johnsonshap <- johnsonshap
    
    assign(id, x, parent.frame())
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
