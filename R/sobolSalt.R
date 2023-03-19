# library(boot)

# LAST MODIFIED: 
# December 9, 2016

# AUTHOR:
# Gilquin Laurent

# SUMMARY:
# Monte Carlo Estimation of Sobol' Indices based on Saltelli schemes.


sobolSalt <- function(model = NULL, X1, X2, scheme="A", nboot = 0, conf = 0.95, ...) {
  
  if ((ncol(X1) != ncol(X2)) | (nrow(X1) != nrow(X2))) {
    stop("The samples X1 and X2 must have the same dimensions")
  }  
  
  if (is.data.frame(X1) & is.data.frame(X2)){
    X1 <- as.matrix(unname(X1))
    X2 <- as.matrix(unname(X2))
  }
  
  else if(!(is.matrix(X1) & is.matrix(X2))){
    stop("The samples X1 and X2 must be matrices or data frames")
  }
  
  if (scheme!="A" & scheme!="B") {
    stop("invalid argument 'scheme', waiting for character 'A' or 'B'")
  }
  
  #creation of the DoE  
  p <- ncol(X1)
  s <- seq(1,p)
  I1 <- (p+1)*s
  
  if(scheme=="A"){
    X0 <- matrix(c(X1,rep(X2,p+1)),ncol=(p+2)*p)
    X0[,I1]=X1[,s]
    X <- matrix(X0[,c(outer(seq(1,(p+2)*p,by=p),seq(0,p-1),'+'))],ncol=p)
  }
  
  if(scheme=="B"){
    X0 <- matrix(c(X1,rep(X2,2*p+1)),ncol=(2*p+2)*p)
    index <- matrix(0,ncol=p-1,nrow=p)
    index[lower.tri(index,diag=F)]=-1
    index <- t(t(index)+s[-1])
    
    I2 <- seq((p+1)*p,2*p^2,by=p)+index
    I <- c(I1,t(I2))
    X0[,I]=X1[,c(s,t(index))]
    X <- matrix(X0[,c(outer(seq(1,(2*p+2)*p,by=p),seq(0,p-1),'+'))],ncol=p)
  }
  
  #creation of the class
  x <- list(model = model, X1 = X1, X2 = X2, scheme=scheme, nboot = nboot,
            conf = conf, X = X, call = match.call()) 
  class(x) <- "sobolSalt"
  
  #calcul of the response for explicit model
  if (! is.null(x$model)) {
    response(x, ...)
    tell(x, ...)
  }  
  return(x)
}


estim.sobolSalt <- function(data, i=1:nrow(data), scheme) {
  
  if(scheme=="A"){
    #quantity of interest
    n <- nrow(data)
    p <- ncol(data)-2
    Y <- matrix(data[i,],nrow=n)
    exclu <- c(1,p+2)
    Mean <- colMeans(Y)
    MeanY <- Mean[1]
    MeanY2 <- (MeanY+Mean[-exclu])/2
    MeanYt <- Mean[p+2]
    MeanY2t <- (MeanYt+Mean[-exclu])/2
    
    a <- rep(0,p)
    b <- rep(0,p)
    c <- rep(0,p)
    out <- rep(0,p)
    
    #first order Sobol' indices
    ind <- rep(1,p)
    ind2 <- 2:(p+1)
    S <- .C("LG_estimator",as.double(Y),as.double(MeanY2),as.integer(p),as.integer(n),as.integer(ind),as.integer(ind2),as.double(a),as.double(b),as.double(c),as.double(out))[[10]]
    #S <- colSums(outer(Y[,1],MeanY2,'-')*t(t(Y[,-exclu])-MeanY2))/(colSums(t(t((Y[,1]+Y[,-exclu])/sqrt(2))-(sqrt(2)+1)*MeanY2)^2)-colSums(Y[,1]*Y[,-exclu]))
    
    #total effect Sobol' indices
    ind <- rep(p+2,p)
    ind2 <- 2:(p+1)
    ST <- 1-.C("LG_estimator",as.double(Y),as.double(MeanY2t),as.integer(p),as.integer(n),as.integer(ind),as.integer(ind2),as.double(a),as.double(b),as.double(c),as.double(out))[[10]]
    #ST <- 1-colSums(outer(Y[,p+2],MeanY2t,'-')*t(t(Y[,-exclu])-MeanY2t))/(colSums(t(t((Y[,p+2]+Y[,-exclu])/sqrt(2))-(sqrt(2)+1)*MeanY2t)^2)-colSums(Y[,p+2]*Y[,-exclu]))

    #storing results
    res <- c(S,ST)
  }
  

  if(scheme=="B"){
    #quantity of interest
    n <- nrow(data)
    p <- (ncol(data)-2)/2
    Y <- matrix(data[i,],nrow=n)
    Mean <- colMeans(Y)
    MeanY <- Mean[1]
    MeanY2 <- (MeanY+Mean[2:(p+1)])/2
    MeanYt <- Mean[2*p+2]
    MeanY2t <- (MeanYt+Mean[2:(p+1)])/2
    MeanY3 <- (MeanYt+Mean[(p+2):(2*p+1)])/2
    MeanY3t <- (MeanY+Mean[(p+2):(2*p+1)])/2
    
    a <- rep(0,p)
    b <- rep(0,p)
    c <- rep(0,p)
    out <- rep(0,p)
    
    #first order Sobol' indices (makes use of double estimates)
    ind <- rep(1,p)
    ind2 <- 2:(p+1)
    S_i <- .C("LG_estimator",as.double(Y),as.double(MeanY2),as.integer(p),as.integer(n),as.integer(ind),as.integer(ind2),as.double(a),as.double(b),as.double(c),as.double(out))[[10]]
    #S_i <- colSums(outer(Y[,1],MeanY2,'-')*t(t(Y[,2:(p+1)])-MeanY2))/(colSums(t(t((Y[,1]+Y[,2:(p+1)])/sqrt(2))-(sqrt(2)+1)*MeanY2)^2)-colSums(Y[,1]*Y[,2:(p+1)]))
    ind <- rep(2*p+2,p)
    ind2 <- (p+2):(2*p+1)
    S_ii <- .C("LG_estimator",as.double(Y),as.double(MeanY3),as.integer(p),as.integer(n),as.integer(ind),as.integer(ind2),as.double(a),as.double(b),as.double(c),as.double(out))[[10]]
    #S_ii <- colSums(outer(Y[,2*p+2],MeanY3,'-')*t(t(Y[,(p+2):(2*p+1)])-MeanY3))/(colSums(t(t((Y[,2*p+2]+Y[,(p+2):(2*p+1)])/sqrt(2))-(sqrt(2)+1)*MeanY3)^2)-colSums(Y[,2*p+2]*Y[,(p+2):(2*p+1)]))
    S <- (S_i+S_ii)/2

    #total effect Sobol' indices (makes use of double estimates)
    ind <- rep(2*p+2,p)
    ind2 <- 2:(p+1)
    ST_i <- .C("LG_estimator",as.double(Y),as.double(MeanY2t),as.integer(p),as.integer(n),as.integer(ind),as.integer(ind2),as.double(a),as.double(b),as.double(c),as.double(out))[[10]]
    #ST_i <- colSums(outer(Y[,2*p+2],MeanY2t,'-')*t(t(Y[,2:(p+1)])-MeanY2t))/(colSums(t(t((Y[,2*p+2]+Y[,2:(p+1)])/sqrt(2))-(sqrt(2)+1)*MeanY2t)^2)-colSums(Y[,2*p+2]*Y[,2:(p+1)]))
    ind <- rep(1,p)
    ind2 <- (p+2):(2*p+1)
    ST_ii <- .C("LG_estimator",as.double(Y),as.double(MeanY3t),as.integer(p),as.integer(n),as.integer(ind),as.integer(ind2),as.double(a),as.double(b),as.double(c),as.double(out))[[10]]
    #ST_ii <- colSums(outer(Y[,1],MeanY3t,'-')*t(t(Y[,(p+2):(2*p+1)])-MeanY3t))/(colSums(t(t((Y[,1]+Y[,(p+2):(2*p+1)])/sqrt(2))-(sqrt(2)+1)*MeanY3t)^2)-colSums(Y[,1]*Y[,(p+2):(2*p+1)]))
    ST <- 1-(ST_i+ST_ii)/2
    
    #second-order Sobol' indices (makes use of double estimates)
    I <- t(combn(p,2))
    MeanY4 <- (Mean[1+I[,1]]+Mean[p+1+I[,2]])/2
    MeanY4b <- (Mean[1+I[,2]]+Mean[p+1+I[,1]])/2
    
    a <- rep(0,nrow(I))
    b <- rep(0,nrow(I))
    c <- rep(0,nrow(I))
    out <- rep(0,nrow(I))
    
    ind <- 1+I[,1]
    ind2 <- p+1+I[,2]
    S2_i <- .C("LG_estimator",as.double(Y),as.double(MeanY4),as.integer(nrow(I)),as.integer(n),as.integer(ind),as.integer(ind2),as.double(a),as.double(b),as.double(c),as.double(out))[[10]]
    #S2_i <- colSums(t(t(Y[,1+I[,1]])-MeanY4)*t(t(Y[,p+1+I[,2]])-MeanY4))/(colSums(t(t((Y[,1+I[,1]]+Y[,p+1+I[,2]])/sqrt(2))-(sqrt(2)+1)*MeanY4)^2)-colSums(Y[,1+I[,1]]*Y[,p+1+I[,2]]))
    ind <- 1+I[,2]
    ind2 <- p+1+I[,1]
    S2_ii <- .C("LG_estimator",as.double(Y),as.double(MeanY4b),as.integer(nrow(I)),as.integer(n),as.integer(ind),as.integer(ind2),as.double(a),as.double(b),as.double(c),as.double(out))[[10]]
    #S2_ii <- colSums(t(t(Y[,1+I[,2]])-MeanY4b)*t(t(Y[,p+1+I[,1]])-MeanY4b))/(colSums(t(t((Y[,1+I[,2]]+Y[,p+1+I[,1]])/sqrt(2))-(sqrt(2)+1)*MeanY4b)^2)-colSums(Y[,1+I[,2]]*Y[,p+1+I[,1]]))
    #closed
    S2 <- (S2_i+S2_ii)/2
    #non-closed
    S2 <- S2-S[I[,1]]-S[I[,2]]
    
    #storing results
    res <- c(S,ST,S2)
  }
  return(res)
}


tell.sobolSalt <- function(x, y = NULL, ...) {
  
  id <- deparse(substitute(x))
  if (! is.null(y)) {
    x$y <- y
  } 
  else if (is.null(x$y)) {
    stop("y not found")
  }
  
  n <- nrow(x$X1)
  p <- ncol(x$X1)
  scheme <- x$scheme
  data <- matrix(x$y, nrow = n)
  
  # output purpose
  rownames <- paste("X",1:p,sep="")
  rownames2 <- paste("X",apply(t(combn(p,2)),1,paste,collapse="X"),sep= "")
  
  x$V <- (n-1)/n*var(data[,1])
  
  # indices estimation + bootstrap
  if (x$nboot == 0) {
    indices <- data.frame("original"=estim.sobolSalt(data, 1:n, scheme))
    x$S <- data.frame("original"=indices[1:p,])
    x$T <- data.frame("original"=indices[(p+1):(2*p),])
    if(scheme=="B"){
      x$S2 <- data.frame("original"=indices[(2*p+1):nrow(indices),])
      rownames(x$S2) <- rownames2
    }
  } else {
    S.boot <- boot(data, estim.sobolSalt, scheme=scheme, R = x$nboot)
    indices <- bootstats(S.boot, x$conf, "basic")
    x$S <- indices[1:p,]
    x$T <- indices[(p+1):(2*p),]
    if(scheme=="B"){
      x$S2 <- indices[(2*p+1):nrow(indices),]
      rownames(x$S2) <- rownames2
    }
  }
  rownames(x$S) <- rownames
  rownames(x$T) <- rownames
  
  assign(id, x, parent.frame())
  return(x)
}


print.sobolSalt <- function(x, ...) {
  
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (! is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    cat("\nModel variance:", x$V, "\n")
    cat("\nFirst order indices:\n")
    print(x$S)
    cat("\nTotal indices:\n")
    print(x$T)
    if(x$scheme=="B"){
      cat("\nSecond order subset indices:\n")
      print(x$S2)
    }
    }
 else {
    cat("(empty)\n")
  }
}

plot.sobolSalt <- function(x, ylim = c(0, 1), choice, ...) {
  
  if (! is.null(x$y)) {
    p <- ncol(x$X)
    pch = c(21, 24)
    if(choice==1){
      nodeplot(x$S, xlim = c(1, p + 1), ylim = ylim, pch = pch[1])
      nodeplot(x$T, xlim = c(1, p + 1), ylim = ylim, labels = FALSE,
               pch = pch[2], at = (1:p)+.3, add = TRUE)
      legend(x = "topright", legend = c("first-order", "total effect"), pch = pch)
    }
    if(x$scheme=="B" & choice==2){
      nodeplot(x$S2, ylim = ylim)
      legend(x = "topright", legend = c("second-order indices"))   
    }
  }
}

ggplot.sobolSalt <- function(data, mapping = aes(), ylim = c(0,1), choice, ..., environment = parent.frame()) {
  x <- data
  
  if (! is.null(x$y)) {
    p <- ncol(x$X)
    pch = c(21, 24)
    if(choice==1){
      g <- nodeggplot(listx = list(x$S,x$T), xname = c("First-order", "Total effect"), ylim = ylim, pch = pch)
    }
    if(x$scheme=="B" & choice==2){
      g <- nodeggplot(listx = list(x$S2), xname = "Second-order indices", ylim = ylim)
    }
    return(g)
  }
}

