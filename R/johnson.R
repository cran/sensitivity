# Johnson indices
#
# Bertrand Iooss and Laura Clouvel 2023

estim.johnson <- function(data, logistic, i = 1:nrow(data)){
  d <- data[i, ]
  
  if (!logistic){
    
    # Computation of the weight matrix
    cor_matrix <- cor(d, use = "pairwise.complete.obs")
    corXX <- cor_matrix[2:ncol(d), 2:ncol(d)] # Correlations between inputs
    corXX.eigen <- eigen(corXX) # Eigenvalues and eigenvectors
    W <- corXX.eigen$vec %*% sqrt(diag(corXX.eigen$val)) %*% t(corXX.eigen$vec)
    
    # Computation of indices (alpha) between output and orthogonal variables
    corXY <- cor_matrix[2:ncol(d), 1] # Correlations between output and inputs
    alpha <- solve(W) %*% corXY # Solve numeric matrix containing coefficients of equation (Ax=B)
    
### With the orthogonalized inputs Z
#     PDQ <- svd(data[,2:ncol(d)])
#      P <- PDQ$u
#      Q <- PDQ$v
#      Z <- P %*% t(Q)
#      alpha <- src(Z,d$Y)$SRC$original
    
    W^2 %*% alpha ^ 2
    
  } else{
    # Tranformation of datas by logistic regression
    xs <- rapply(d[, 2:ncol(d)],scale,c("numeric"),how="replace")
    xs <- data.frame(xs)
    m <- glm(y~., data=data.frame(y = d[,1], xs), family="binomial", maxit=100)
    R0 <-1-m$deviance/m$null.deviance
    gp_estimated <- data.frame(m$linear.predictors)
    gp_scale <- rapply(gp_estimated,scale,c("numeric"),how="replace")
    dlogit <- data.frame(Y = gp_scale, xs)
    
    # Computation of the weight matrix
    cor_matrix <- cor(dlogit, use = "pairwise.complete.obs")
    corXX <- cor_matrix[2:ncol(dlogit), 2:ncol(dlogit)] # Correlations between inputs
    corXX.eigen <- eigen(corXX) # Eigenvalues and eigenvectors
    W <- corXX.eigen$vec %*% sqrt(diag(corXX.eigen$val)) %*% t(corXX.eigen$vec)
    
    # Computation of indices (alpha) between output and orthogonal variables
    corXY <- cor_matrix[2:ncol(dlogit), 1] # Correlations between output and inputs
    alpha <- solve(W) %*% corXY # Solve numeric matrix containing coefficients of equation (Ax=B)
    
    W^2 %*% alpha ^ 2 * R0
  }
}


johnson <- function(X, y, rank = FALSE, logistic = FALSE, nboot = 0, conf = 0.95) {
  data <- data.frame(Y = y, X)
  
  if (logistic) rank <- FALSE # Impossible to perform logistic regression with a rank transformation
  
  if (rank) {
    for (i in 1:ncol(data)) {
      data[,i] <- rank(data[,i])
    }
  }
  
  if (nboot == 0) {
    johnson <- data.frame(original = estim.johnson(data, logistic ))
    rownames(johnson) <- colnames(X)
  } else {
    boot.johnson <- boot(data, estim.johnson, logistic = logistic, R = nboot)
    johnson <- bootstats(boot.johnson, conf, "basic")
    rownames(johnson) <- colnames(X)
  }

  out <- list(X = X, y = y, rank = rank, nboot = nboot, conf = conf,
              call = match.call())
  class(out) <- "johnson"
  out$johnson <- johnson
    
  return(out)
}


print.johnson <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  cat("\nJohnson indices:\n")
  print(x$johnson)
}


plot.johnson<- function(x, ylim = c(0,1), ...) {  
  nodeplot(x$johnson, ylim = ylim, main = "Johnson indices")
}

ggplot.johnson <- function(data, mapping = aes(), ylim = c(0,1), ..., environment = parent.frame()) {  
  x <- data
  nodeggplot(listx = list(x$johnson), xname = "Johnson indices", ylim = ylim, title = "Johnson indices")
}
