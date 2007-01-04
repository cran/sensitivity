
## File: sobol_saltelli02.R
## Description: The method of Sobol, estimation sheme from Saltelli (2002)
## Author: Gilles Pujol


sobol.saltelli02 <- function(model = NULL, x1, x2, nboot = 0, conf = 0.95, ...){
  # ARGUMENTS

  if ( (ncol(x1) != ncol(x2)) | (nrow(x1) != nrow(x2)) )
    stop("The two samples x1 and x2 must have the same dimensions")
  p <- ncol(x1)
  
  # DESIGN OF EXPERIMENTS
  
  x <- rbind(x1, x2)
  for ( i in 1:p ){
    xb <- x1
    xb[, i] <- x2[, i]
    x <- rbind(x, xb) 
  }

  # OBJECT OF CLASS "sobol.saltelli02"
  
  sa <- list(model = model, x1 = x1, x2 = x2, nboot = nboot, conf = conf, x = x,
             y = NULL, S1 = NULL, St = NULL, call = match.call())

  class(sa) <- "sobol.saltelli02"

  # COMPUTATION OF THE RESPONSE AND SENSITIVITY ANALYSIS
  
  if (!is.null(sa$model)){
      response(sa, ...)
      tell(sa)
  }

  # RETURN OF THE OBJECT OF CLASS "sobol"
  
  return(sa)
}


estim.sobol.saltelli02 <- function(data, i = 1 : nrow(data)){
  d <- as.matrix(data[i, ]) # as.matrix pour que colSums renvoie un numeric...
  n <- nrow(d)
  
  V <- var(d[, 1])
  first.order <- (colSums(d[, - c(1, 2)] * d[, 2]) / (n - 1) -
                  mean(d[,1] * d[, 2])) / V
  total.order <- 1 - (colSums(d[, - c(1, 2)] * d[, 1]) / (n - 1) -
                  mean(d[, 1])^2) / V
  
  return(c(first.order, total.order))
}


tell.sobol.saltelli02 <- function(sa, y = NULL){
  id <- deparse(substitute(sa))
  n <- nrow(sa$x1)
  p <- ncol(sa$x1)
  
  # EXTERNAL MODEL

  if (! is.null(y))
    sa$y <- y
  
  # ESTIMATION OF THE INDICES

  data <- matrix(sa$y, nr = n, nc = p + 2)

  if (sa$nboot == 0){
    # single estimation
    S <- data.frame(original = estim.sobol.saltelli02(data))
  }
  else{
    # bootstrap estimation
    boot.sobol <- boot(data, estim.sobol.saltelli02, R = sa$nboot)
    S <- bootstats(boot.sobol, sa$conf, "basic")
  }

  sa$S1 <- subset(S, subset = c(rep(TRUE, p), rep(FALSE, p)))
  rownames(sa$S1) <- colnames(sa$x1)

  sa$St <- subset(S, subset = c(rep(FALSE, p), rep(TRUE, p)))
  rownames(sa$St) <- colnames(sa$x1)
    
  # SAVING OF THE INDICES
  
  assign(id, sa, parent.frame())
}


print.sobol.saltelli02 <- function(x, ...){
  cat("\nSOBOL METHOD\n")
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  cat("\nModel runs:", length(x$y), "\n")
  cat("\nFirst order indices:\n")
  print(x$S1)
  cat("\nTotal indices:\n")
  print(x$St)
}


plot.sobol.saltelli02 <- function(x, ...){
  pch <- c(21, 24)
  p <- ncol(x$x1)
  nodeplot(x$S1, xlim = c(1,p+1), ylim = c(0,1), labels = colnames(x$x1),
           pch = pch[1])
  nodeplot(x$St, xlim = c(1,p+1), ylim = c(0,1), labels = FALSE,
           pch = pch[2], at = (1:p)+.3, add = TRUE)
  legend(x = "topright", legend = c("first order", "total order"), pch = pch)
}
