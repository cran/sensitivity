# Sobol' indices estimation (Saltelli 2002)
#
# Gilles Pujol 2006


sobol2002 <- function(model = NULL, X1, X2, nboot = 0, conf = 0.95, ...) {
  if ((ncol(X1) != ncol(X2)) | (nrow(X1) != nrow(X2)))
    stop("The samples X1 and X2 must have the same dimensions")
  p <- ncol(X1)
  
  X <- rbind(X1, X2)
  for (i in 1:p) {
    Xb <- X1
    Xb[,i] <- X2[,i]
    X <- rbind(X, Xb) 
  }
  
  x <- list(model = model, X1 = X1, X2 = X2, nboot = nboot, conf = conf, X = X,
            call = match.call())
  class(x) <- "sobol2002"
  
  if (!is.null(x$model)) {
    response(x, ...)
    tell(x)
  }
  
  return(x)
}


estim.sobol2002 <- function(data, i = 1 : nrow(data)) {
  d <- as.matrix(data[i, ]) # as.matrix for colSums
  n <- nrow(d)
  V <- var(d[, 1])
  VCE <- (colSums(d[, - c(1, 2)] * d[, 2]) / (n - 1) - mean(d[,1] * d[,2]))
  VCE.compl <- (colSums(d[, - c(1, 2)] * d[, 1]) / (n - 1) - mean(d[, 1])^2)
  c(V, VCE, VCE.compl)
}


tell.sobol2002 <- function(x, y = NULL, return.var = NULL, ...) {
  id <- deparse(substitute(x))

  if (! is.null(y)) {
    x$y <- y
  } else if (is.null(x$y)) {
    stop("y not found")
  }

  p <- ncol(x$X1)
  n <- nrow(x$X1)

  data <- matrix(x$y, nrow = n)

  # estimation of the partial variances (V, D1 and Dt)
  
  if (x$nboot == 0){
    V <- data.frame(original = estim.sobol2002(data))
  }
  else{
    V.boot <- boot(data, estim.sobol2002, R = x$nboot)
    V <- bootstats(V.boot, x$conf, "basic")
  }
  rownames(V) <- c("global", colnames(x$X1), paste("-", colnames(x$X1), sep = ""))

  # estimation of the Sobol' indices (S1 and St)

  if (x$nboot == 0) {
    S <- V[2:(p + 1), 1, drop = FALSE] / V[1,1]
    T <- 1 - V[(p + 2):(2 * p + 1), 1, drop = FALSE] / V[1,1]
  } else {
    S.boot <- V.boot
    S.boot$t0 <- V.boot$t0[2:(p + 1)] / V.boot$t0[1]
    S.boot$t <- V.boot$t[,2:(p + 1)] / V.boot$t[,1]
    S <- bootstats(S.boot, x$conf, "basic")
    
    T.boot <- V.boot
    T.boot$t0 <- 1 - V.boot$t0[(p + 2):(2 * p + 1)] / V.boot$t0[1]
    T.boot$t <- 1 - V.boot$t[,(p + 2):(2 * p + 1)] / V.boot$t[,1]
    T <- bootstats(T.boot, x$conf, "basic")
  }
  rownames(S) <- colnames(x$X1)
  rownames(T) <- colnames(x$X1)

  # return
  x$V <- V
  x$S <- S
  x$T <- T

  for (i in return.var) {
    x[[i]] <- get(i)
  }

  assign(id, x, parent.frame())
}


print.sobol2002 <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (!is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    cat("\nFirst order indices:\n")
    print(x$S)
    cat("\nTotal indices:\n")
    print(x$T)
  }
}


plot.sobol2002 <- function(x, ylim = c(0, 1), ...) {
  if (!is.null(x$y)) {
    p <- ncol(x$X1)
    pch = c(21, 24)
    nodeplot(x$S, xlim = c(1, p + 1), ylim = ylim, pch = pch[1])
    nodeplot(x$T, xlim = c(1, p + 1), ylim = ylim, labels = FALSE,
             pch = pch[2], at = (1:p)+.3, add = TRUE)
    legend(x = "topright", legend = c("main effect", "total effect"), pch = pch)
  }
}

ggplot.sobol2002 <- function(data, mapping = aes(), ylim = c(0,1), ..., environment = parent.frame()) {
  x <- data
  if (!is.null(x$y)) {
    p <- ncol(x$X1)
    pch = c(21, 24)
    nodeggplot(listx = list(x$S,x$T), xname = c("Main effet","Total effect"), ylim = ylim, pch = pch)
  }
}

plotMultOut.sobol2002 <- function(x, ylim = c(0, 1), ...) {
  if (!is.null(x$y)) {
    p <- ncol(x$X1)
    if (!x$ubiquitous){
      stop("Cannot plot functional indices since ubiquitous option was not activated")
    }else{
      if (x$Tot == T) par(mfrow=c(2,1))
      plot(0,ylim=ylim,xlim=c(1,x$q),main="First order Sobol indices",ylab="",xlab="",type="n")
      for (i in 1:p) lines(x$Sfct[,i],col=i)
      legend(x = "topright", legend = dimnames(x$X1)[[2]], lty=1, col=1:p, cex=0.6)
      
      if (x$Tot == T){
        plot(0,ylim=ylim,xlim=c(1,x$q),main="Total Sobol indices",ylab="",xlab="",type="n")
        for (i in 1:p) lines(x$Tfct[,i],col=i)
        legend(x = "topright", legend = dimnames(x$X1)[[2]], lty=1, col=1:p, cex=0.6)
      }
    }
  }
}
