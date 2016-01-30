# Sobol' indices estimation (Martinez 2011)
# Plus: Theoretical confidence intervals from correlation coefficient-based confidence interval
#
# J-M. Martinez, Analyse de sensibilite globale par decomposition de la variance,
# Presentation a la journee des GdR Ondes et MASCOT-NUM, 13 janvier 2011, 
# Institut Henri Poincare, Paris, France.
#
# Bertrand Iooss (2015)


sobolmartinez <- function(model = NULL, X1, X2, nboot = 0, conf = 0.95, ...) {
  if ((ncol(X1) != ncol(X2)) | (nrow(X1) != nrow(X2)))
    stop("The samples X1 and X2 must have the same dimensions")
  p <- ncol(X1)
  
  X <- rbind(X1,X2)
    for (i in 1:p) {
    Xb <- X1
    Xb[,i] <- X2[,i]
    X <- rbind(X, Xb) 
  }
  
  x <- list(model = model, X1 = X1, X2 = X2, nboot = nboot, conf = conf, X = X,
            call = match.call())
  class(x) <- "sobolmartinez"
  
  if (!is.null(x$model)) {
    response(x, ...)
    tell(x)
  }
  
  return(x)
}


estim.sobolmartinez <- function(data, i = 1 : nrow(data), estimStd=FALSE, conf=0) {
  d <- as.matrix(data[i, ]) # as.matrix for colSums
  n <- nrow(d)
  p <- ncol(d)-2
  
  V <- var(d[, 1])
  ecor <- rep(0,p) ; ecorcompl <- rep(0,p)
  if(estimStd) {
    VV <- matrix(V,nrow=1,ncol=3,dimnames=list(1,c("estim","CIinf","CIsup")))
    estcor <- matrix(0,nrow=p,ncol=3,dimnames=list(2:(p+1),c("estim","CIinf","CIsup")))
    estcorcompl <- matrix(0,nrow=p,ncol=3,dimnames=list((p+2):(2*p+1),c("estim","CIinf","CIsup")))
  }
  for(ii in 1:p) {
      ecor[ii] <- cor(d[,2],d[,ii+2],use="pairwise.complete.obs")
      ecorcompl[ii] <- cor(d[,1],d[,ii+2],use="pairwise.complete.obs")
    
    if(estimStd) {
      estcor[ii,1] <- ecor[ii]
      confcor <- cor.test(d[,2],d[,ii+2],conf.level=conf)
      estcor[ii,2] <- confcor$conf.int[1]
      estcor[ii,3] <- confcor$conf.int[2]
      estcorcompl[ii,1] <- ecorcompl[ii]
      confcor <- cor.test(d[,1],d[,ii+2],conf.level=conf)
      estcorcompl[ii,2] <- confcor$conf.int[2]
      estcorcompl[ii,3] <- confcor$conf.int[1] # on intervertit car apres on prend l'oppose
    }
  }
  if(estimStd) 
    { return(rbind(VV, estcor, estcorcompl))}
  else 
    { return(c(V, ecor, ecorcompl))}
}


tell.sobolmartinez <- function(x, y = NULL, return.var = NULL, ...) {
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
    V <- data.frame(original = estim.sobolmartinez(data, 1:n, TRUE, x$conf))
    colnames(V) <- c("original","min. c.i.","max. c.i.")
  }
  else{
    V.boot <- boot(data, estim.sobolmartinez, R = x$nboot)
    V <- bootstats(V.boot, x$conf, "basic")
    rownames(V) <- c("global", colnames(x$X1), paste("-", colnames(x$X1), sep = ""))
  }
  

  # estimation of the Sobol' indices (S1 and St)

  if (x$nboot == 0) {
   S <- V[2:(p + 1), 1:3, drop = FALSE]
   T <- 1 - V[(p + 2):(2 * p + 1), 1:3, drop = FALSE]
   
  } else {
    S.boot <- V.boot
    S.boot$t0 <- V.boot$t0[2:(p + 1)]
    S.boot$t <- V.boot$t[,2:(p + 1)]
    S <- bootstats(S.boot, x$conf, "basic")
    
    T.boot <- V.boot
    T.boot$t0 <- 1 - V.boot$t0[(p + 2):(2 * p + 1)]
    T.boot$t <- 1 - V.boot$t[,(p + 2):(2 * p + 1)]
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


print.sobolmartinez <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (!is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    cat("\nFirst order indices:\n")
    if(!is.null(x$S$original.CIinf)) {
      df=data.frame(pointEstimate=x$S[,1],  minCI=x$S[,2], maxCI=x$S[,3])
      colnames(df)=c("estimate","min. c.i.","max. c.i.")
      print(df)
    } else {
      print(x$S)
    }
    cat("\nTotal indices:\n")
    if(!is.null(x$T$original.CIinf)) {
      df=data.frame(pointEstimate=x$T[,1], minCI=x$T[,2], maxCI=x$T[,3])
      colnames(df)=c("estimate","min. c.i.","max. c.i.")
      print(df)
    } else {
      print(x$T)
    }
  }
}


plot.sobolmartinez <- function(x, ylim = c(0, 1), ...) {
  if (!is.null(x$y)) {
    p <- ncol(x$X1)
    pch = c(21, 24)
    nodeplot(x$S, xlim = c(1, p + 1), ylim = ylim, pch = pch[1])
    nodeplot(x$T, xlim = c(1, p + 1), ylim = ylim, labels = FALSE,
             pch = pch[2], at = (1:p)+.3, add = TRUE)
    legend(x = "topright", legend = c("main effect", "total effect"), pch = pch)
  }
}