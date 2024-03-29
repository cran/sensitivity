# Asymptotically efficient Sobol' indices estimation (Monod 2006, Janon et al. 2014)
#
# Alexandre Janon 2012
# Laurent Gilquin 2014


sobolEff <- function(model = NULL, X1, X2, order=1, nboot = 0, conf = 0.95, ...) {
  
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
  
  if(!(order!=0 | order!=1 | order!=2)){
    stop("invalid argument 'order', waiting for value 0, 1 or 2.")
  }
  
  #creation of the DoE
  
  #cas order 1
  if(order==1){
    p <- ncol(X1)
    X0 <- matrix(c(X1,rep(X2,p)),ncol=(p+1)*p)
    X0[,(p+1)*seq(1,p)]=X1[,seq(1,p)]
    X <- matrix(X0[,c(outer(seq(1,(p+1)*p,by=p),seq(0,p-1),'+'))],ncol=p)
  }
  
  #cas total
  if(order==0){
    p <- ncol(X1)
    X0 <- matrix(rep(X1,p+1),ncol=(p+1)*p)
    X0[,(p+1)*seq(1,p)]=X2[,seq(1,p)]
    X <- matrix(X0[,c(outer(seq(1,(p+1)*p,by=p),seq(0,p-1),'+'))],ncol=p)
  }
  
  #cas order 2
  if(order==2){
    p <- ncol(X1) 
    i <- 1
    j <- 1
    I <- matrix(0,ncol=2,nrow=p*(p-1)/2)
    for (l in 1:(p*(p-1)/2)) {
      i <- i+(j==p)
      j <- ((j==p)*i + (j<p)*j)+1
      I[l,]=c(i,j)
    }
    X0 <- matrix(c(X1,rep(X2,p*(p-1)/2)),ncol=p*(1+p*(p-1)/2))
    X0[,c(t(seq(p,p^2*(p-1)/2,by=p)+I))]=X1[,c(t(I))]
    X <- matrix(X0[,c(outer(seq(1,p*(1+p*(p-1)/2),by=p),seq(0,p-1),'+'))],ncol=p)
  }
  
  #creation of the class
  x <- list(model = model, X1 = X1, X2 = X2, order = order, nboot = nboot,
            conf = conf, X = X, call = match.call()) 
  class(x) <- "sobolEff"
  
  #calcul of the response for explicit model
  if (! is.null(x$model)) {
    response(x, ...)
    x=tell(x, ...)
  }  
  return(x)
}


estim.sobolEff <- function(data, i=1:nrow(data), estimStd, conf) {
  
  #quantity of interest
  n <- nrow(data)
  Y <- as.matrix(data[i,],nrow=n)
  Mean <- colMeans(Y)
  MeanY <- Mean[1]
  MeanY2 <- (MeanY+Mean[-1])/2
  
  p <- ncol(data)
  a <- rep(0,p-1)
  b <- rep(0,p-1)
  c <- rep(0,p-1)
  out <- rep(0,p-1)
  
  #Sobol' indices
  ind <- rep(1,p-1)
  ind2 <- 2:p
  S <- .C("LG_estimator",as.double(Y),as.double(MeanY2),as.integer(p-1),as.integer(n),as.integer(ind),as.integer(ind2),as.double(a),as.double(b),as.double(c),as.double(out))[[10]]
  #S <- colSums(outer(Y[,1],MeanY2,'-')*(t(t(Y[,-1])-MeanY2)))/(colSums(t(t((Y[,1]+Y[,-1])/sqrt(2))-(sqrt(2)+1)*MeanY2)^2)-colSums(Y[,1]*Y[,-1]))
  res <- S
  
  #estimation of the standard deviation of the Sobol' indice estimators
  if(estimStd){
    term1 <- Y[,1]-MeanY
    term2 <- Y[,-1]-MeanY
    U <- term1*term2-1/2*t(S*t(term1^2+term2^2))
    std <- (n/(n-1))*apply(U,2,var)/(var(Y[,1])^2)
    stderr <- sqrt(std/n)
    R <- qnorm((1+conf)/2)*stderr
    
    res <- list(S=round(S,6),stderr=round(stderr,6),Rad=round(R,6))
  }
  return(res)
}


tell.sobolEff <- function(x, y = NULL, ...) {
  
  id <- deparse(substitute(x))
  if (! is.null(y)) {
    x$y <- y
  } 
  else if (is.null(x$y)) {
    stop("y not found")
  }
  
  n <- nrow(x$X1)
  data <- matrix(x$y, nrow = n)
  x$V <- (n-1)/n*var(data[,1])
  
  if (x$nboot == 0) {
    res <- estim.sobolEff(data, 1:n, estimStd=TRUE, x$conf)
    x$S <- data.frame(res$S,res$stderr,res$S-res$Rad,res$S+res$Rad)
    colnames(x$S)=c("original", "std. error" , "min. c.i." , "max. c.i.")
  } 
  else {
    S.boot <- boot(data, estim.sobolEff, estimStd=FALSE, R = x$nboot)
    x$S <- bootstats(S.boot, x$conf, "basic")
  }
  if(x$order==0){
    x$S[c("original","min. c.i.","max. c.i.")] <- 1-x$S[c("original","max. c.i.","min. c.i.")]
  }
  
  d <- ncol(x$X1)
  if (x$order==1 | x$order==0){
    rownames <- paste("X",1:d,sep="")
  }
  if (x$order==2){
    rownames <- paste("X",apply(t(combn(d,2)),1,paste,collapse=""),sep= "")
  }
  rownames(x$S) <- rownames
  
  assign(id, x, parent.frame())
  return(x)
}


print.sobolEff <- function(x, ...) {
  
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (! is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    cat("\nModel variance:", x$V, "\n")
    if (! is.null(x$S)) {
      cat("\n\n\nSobol indices\n")
      print(x$S)
    }
  } else {
    cat("(empty)\n")
  }
}


plot.sobolEff <- function(x, ylim = c(0, 1), ...) {
  
  if (! is.null(x$y)) {
    nodeplot(x$S, ylim = ylim)
    if (x$order==1){
      legend(x = "topright", legend = c("First order indices"))
    }
    if (x$order==2){
      legend(x = "topright", legend = c("Second order subset indices"))
    }
    if (x$order==0){
      legend(x = "topright", legend = c("Total indices"))
    }
  }
}

ggplot.sobolEff <- function(data, mapping = aes(), ylim = c(0,1), ..., environment = parent.frame()) {
  x <- data
  
  if (! is.null(x$y)) {
    if (x$order==1){
      title <- "First order indices"
    }
    if (x$order==2){
      title <- "Second order subset indices"
    }
    if (x$order==0){
      title <- "Total indices"
    }
    nodeggplot(listx = list(x$S), xname = title, title = title, ylim = ylim)
  }
}

