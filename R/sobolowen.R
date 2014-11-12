# Sobol' indices estimation by the method of Owen
#Authors: B. Ramos and T. Touati (2014)

sobolowen <- function(model = NULL, X1, X2, X3, nboot = 0, conf = 0.95, varest = 2, ...) {
  if ((ncol(X1) != ncol(X2)) | (nrow(X1) != nrow(X2)) | (ncol(X2) != ncol(X3)) | (nrow(X2) != nrow(X3)))
    stop("The samples X1, X2 and X3 must have the same dimensions")
  p <- ncol(X1)
  
  X <- rbind(X1, X2)
  for (i in 1:p) {
    Xb <- X1
    Xb[,i] <- X3[,i]
    X <- rbind(X, Xb)
  }
  for (i in 1:p) {
    Xb <- X2
    Xb[,i] <- X1[,i]
    X <- rbind(X, Xb)
  }
  for (i in 1:p) {
    Xb <- X3
    Xb[,i] <- X2[,i]
    X <- rbind(X, Xb)
  }
  
  x <- list(model = model, X1 = X1, X2 = X2, X3 = X3, nboot = nboot, conf = conf, X = X,
            call = match.call())
  class(x) <- "sobolowen"
  
  if (!is.null(x$model)) {
    response(x, ...)
    tell(x,varest=varest)
  }
  
  return(x)
}

response <- function(x, loop = FALSE, ...) {
  id <- deparse(substitute(x))
  if (class(x$model) == "function") {
    if (loop) {
      n <- nrow(x$X)
      y <- numeric(n)
      for (i in 1:n) {
        y[i] <- x$model(x$X[i,], ...)
      }
    } else {
      y <- x$model(x$X, ...)
    }
  } else if (TRUE %in% (paste("predict.", class(x$model), sep="") %in% methods(predict))) {
    y <- predict(x$model, x$X, ...)
  } else {
    stop("The model isn't a function or does not have a predict method")
  }
  
  if (class(y) != "numeric") {
    y <- as.numeric(y)
    warning("Conversion of the response to numeric")
  }
  
  x$y <- y
  assign(id, x, parent.frame())
} 

estim.sobolowen <- function(data, i=1:nrow(data), varest=2) {
  d <- as.matrix(data[i, ]) # as.matrix for colSums
  n <- nrow(d)
  p <- (ncol(d)-2)/3
  if (varest==1) {
    V <- var(d[, 1])    
  } else {
    V <- numeric(0)    
    for (k in 1:p) {
      V[k] <- mean(apply(d[,c(1,2,2+k,2+k+p)]^2,1,mean))-(mean(apply(d[,c(1,2,2+k,2+k+p)],1,mean)))^2
    }
  }
  VCE <- (colSums((d[,1] - d[, 3:(2+p)])*(d[, (3+p):(2+2*p)] - d[,2])) / n)
  VCE.compl <- V - (colSums((d[,2] - d[, (3+2*p):(2+3*p)])*(d[, (3+p):(2+2*p)] - d[,1])) / n)
  # COMM: VCE correspond a underline{tau_u}=\int (f(Y)-f(Y_U:Z_{-U}))(f(X_U:Y_{-U})-f(X))
  # COMM: VCE.compl CORRESPOND A overline{tau_u}=sigma^2-underline{tau_{-u}}=sigma^2 - \int (f(Y)-f(Y_U:Z_{-U}))(f(X_U:Y_{-U})-f(X)) qui est obtenu en echangeant x et y de l'estimateur d'en haut et en remplacant u par -u
  c(V, VCE, VCE.compl)
}

tell.sobolowen <- function(x, y = NULL, return.var = NULL, varest=2, ...) {
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
    V <- data.frame(original = estim.sobolowen(data,varest=varest))
  }
  else{
    #    eval(parse(text=paste("func <- function(data,i) {estim.sobolowen(data,i,varest=",varest,") }",sep="")))
    func <- function(data,i) {estim.sobolowen(data,i,varest)}
    V.boot <- boot(data=data, statistic=func, R = x$nboot)
    V <- bootstats(V.boot, x$conf, "basic")
  }
  rownames(V) <- c(paste("global",colnames(x$X1)), colnames(x$X1), paste("-", colnames(x$X1), sep = ""))
  
  # estimation of the Sobol' indices (S1 and St)
  if (varest==1) {
    k <- 1
  } else {
    k <- p
  }
  if (x$nboot == 0) {
    S <- V[(k + 1):(p + k), 1, drop = FALSE] / V[1:k,1]
    T <- V[(p + k + 1):(2 * p + k), 1, drop = FALSE] / V[1:k,1]
  } else {
    S.boot <- V.boot
    S.boot$t0 <- V.boot$t0[(k+1):(p + k)] / V.boot$t0[1:k]
    S.boot$t <- V.boot$t[,(k+1):(p + k)] / V.boot$t[,(1:k)]
    S <- bootstats(S.boot, x$conf, "basic")
    
    T.boot <- V.boot
    T.boot$t0 <- V.boot$t0[(p + k + 1):(2 * p + k)] / V.boot$t0[1:k]
    T.boot$t <- V.boot$t[,(p + k + 1):(2 * p + k)] / V.boot$t[,(1:k)]
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

print.sobolowen <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (!is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    cat("\nFirst order indices:\n")
    print(x$S)
    cat("\nTotal indices:\n")
    print(x$T)
  }
}

plot.sobolowen <- function(x, ylim = c(0, 1), ...) {
  if (!is.null(x$y)) {
    p <- ncol(x$X1)
    pch = c(21, 24)
    nodeplot(x$S, xlim = c(1, p + 1), ylim = ylim, pch = pch[1])
    nodeplot(x$T, xlim = c(1, p + 1), ylim = ylim, labels = FALSE,
             pch = pch[2], at = (1:p)+.3, add = TRUE)
    legend(x = "topright", legend = c("main effect", "total effect"), pch = pch)
  }
}


nodeplot <- function (x, xlim = NULL, ylim = NULL, labels = TRUE, col = par("col"), 
                      pch = 21, bg = "white", add = FALSE, at = NULL, ...) 
{
  n <- nrow(x)
  if (is.null(xlim)) {
    xlim <- c(1, n)
  }
  if (is.null(ylim)) {
    ylim <- c(min(x), max(x))
  }
  if (is.null(at)) {
    at <- 1:n
  }
  if (add) {
    par(new = TRUE)
  }
  plot(0, xlim = xlim, ylim = ylim, axes = FALSE, xlab = "", 
       ylab = "", type = "n", ...)
  if (class(labels) == "logical") {
    if (labels) {
      axis(side = 1, at = at, labels = rownames(x))
    }
    else {
      axis(side = 1, at = at, labels = FALSE, tick = FALSE)
    }
  }
  else if (class(labels) == "character") {
    axis(side = 1, at = at, labels = labels)
  }
  axis(side = 2)
  box()
  if ("bias" %in% colnames(x)) {
    xx <- x[["original"]] - x[["bias"]]
  }
  else {
    xx <- x[["original"]]
  }
  if (("min. c.i." %in% colnames(x)) & "max. c.i." %in% colnames(x)) {
    for (i in 1:n) {
      lines(c(at[i], at[i]), c(x[["min. c.i."]][i], x[["max. c.i."]][i]), 
            col = col)
    }
  }
  points(at, xx, col = col, pch = pch, bg = bg)
}

#SOBOLOWEN
bootstats <- function (b, conf = 0.95, type = "norm") 
{
  p <- length(b$t0)
  lab <- c("original", "bias", "std. error", "min. c.i.", "max. c.i.")
  out <- as.data.frame(matrix(nrow = p, ncol = length(lab), 
                              dimnames = list(NULL, lab)))
  for (i in 1:p) {
    out[i, "original"] <- b$t0[i]
    out[i, "bias"] <- mean(b$t[, i]) - b$t0[i]
    out[i, "std. error"] <- sd(b$t[, i])
    if (type == "norm") {
      ci <- boot.ci(b, index = i, type = "norm", conf = conf)
      if (!is.null(ci)) {
        out[i, "min. c.i."] <- ci$norm[2]
        out[i, "max. c.i."] <- ci$norm[3]
      }
    }
    else if (type == "basic") {
      ci <- boot.ci(b, index = i, type = "basic", conf = conf)
      if (!is.null(ci)) {
        out[i, "min. c.i."] <- ci$basic[4]
        out[i, "max. c.i."] <- ci$basic[5]
      }
    }
  }
  return(out)
}
