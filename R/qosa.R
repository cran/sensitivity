# Quantile-oriented sensitivity analysis 
# Sebastien Da Veiga (2020)


rkde.w <- function(n, fhat)
{
  x <- fhat$x
  nsamp <- length(x)
  w <- fhat$w
  x.ind <- sample(1:nsamp, size=n, replace=TRUE, prob=w)
  h <- fhat$h
  rkde.val <- x[x.ind] + h*rnorm(n=n, mean=0, sd=1)
  return(rkde.val)
}

qosa <- function(model = NULL, X1, X2 = NULL, type = "quantile", alpha = 0.1, split.sample = 2/3, nsample = 1e4, nboot = 0, conf = 0.95, ...) {
  
  if (is.null(X2)){
    n <- nrow(X1)
    n1 <- floor(split.sample*n)
    X2 <- X1[(n1+1):n,]
    X1 <- X1[1:n1,]
  }
  
  if (ncol(X1) != ncol(X2)) {
    stop("The samples X1 and X2 must have the same dimensions")
  } 
  
  if (is.data.frame(X1) & is.data.frame(X2)){
    X1 <- as.matrix(X1)
    X2 <- as.matrix(X2)
  }else if(!(is.matrix(X1) & is.matrix(X2))){
    stop("The samples X1 and X2 must be matrices or data frames")
  }
  
  if (type!="quantile" & type!="mean"){
    stop("Invalid argument 'type', should be either 'quantile' or 'mean'")
  }
  
  X <- rbind(X1,X2)
  
  x <- list(model = model, X1 = X1, X2 = X2, type = type, alpha = alpha, nsample = nsample, nboot = nboot,
            conf = conf, X = X, call = match.call()) 
  class(x) <- "qosa"
  
  #calcul of the response for explicit model
  if (! is.null(x$model)) {
    response(x, ...)
    x=tell(x, ...)
  }  
  return(x)
}

estim.qosa <- function(data, i=1:nrow(data), type, alpha, n1, nsample) {
  ptot <- ncol(data)
  p <- ptot - 1
  id1 <- i[i<=n1] # indices for first sample
  id2 <- setdiff(i,id1) # indices for second sample
  X1 <- data[id1,1:p]
  Y1 <- data[id1,ptot]
  n1 <- nrow(X1)
  X2 <- data[id2,1:p]
  Y2 <- data[id2,ptot]
  n2 <- nrow(X2)
  S = matrix(0,nrow=p,ncol=1)
  if (type=="quantile"){
    qY <- quantile(Y1,alpha) # compute quantile of output with sample Y1
  }
  EY <- mean(Y2) # compute mean of output with sample Y2
  
  for (j in 1:p){
    # Step 1: estimate diagonal bandwidth  with sample (X1,Y1)
    h <- diag(ks::Hpi.diag(cbind(X1[,j],Y1)))
    # Step 2: loop on samples X2 to compute conditional quantiles or means at each point of X2
    if (type=="quantile"){
    cond.summary <- sapply(1:n2,function(k){
      # Build kde object corresponding to conditional pdf
      # Step 2.1: reconstruct kernel matrix on X
      w.cond <- exp(-0.5*(X2[k,j]-X1[,j])^2/h[1])
      if (max(w.cond)>0){
        w.cond <- w.cond/sum(w.cond)*n1
        # Step 2.2: conditional pdf estimate
        fhat.cond <- ks::kde(Y1,H=h[2],w=w.cond)
        # Step 2.3: compute quantile
        # We first sample from the conditonal pdf 
        # and compute the empirical quantile on this sample
        quantile(rkde.w(nsample,fhat.cond),alpha)
      }else{
        NA
      }
    })
    }else{
      cond.summary <- sapply(1:n2,function(k){
        # Build kde object corresponding to conditional pdf
        # Step 2.1: reconstruct kernel matrix on X
        w.cond <- exp(-0.5*(X2[k,j]-X1[,j])^2/h[1])
        if (max(w.cond)>0){
          w.cond <- w.cond/sum(w.cond)*n1
          # Step 2.2: conditional pdf estimate
          fhat.cond <- ks::kde(Y1,H=h[2],w=w.cond)
          # Step 2.3: compute quantile
          # We first sample from the conditonal pdf 
          # and compute the empirical mean on this sample (equivalent to Nadaraya-Watson estimate)
          mean(rkde.w(nsample,fhat.cond))
        }else{
          NA
        }
      })
    }
    id.notna <- which(!is.na(cond.summary))
    Yfilter <- Y2[id.notna]
    cond.filter <- cond.summary[id.notna]
    n.filter <- length(id.notna)
    if (type=="quantile"){
      R <- Yfilter*(as.numeric(Yfilter>cond.filter)/(1-alpha)-1)
      Z <- Yfilter*(as.numeric(Yfilter>qY)/(1-alpha)-1)
      S[j] <- 1 - mean(R)/mean(Z)
    }else{
      S[j] <- var(cond.filter)/var(Yfilter)
    }
  }
  return(S)
}

tell.qosa <- function(x, y = NULL, ...) {
  
  id <- deparse(substitute(x))
  if (! is.null(y)) {
    x$y <- y
  } 
  else if (is.null(x$y)) {
    stop("y not found")
  }
  
  n <- nrow(x$X)
  n1 <- nrow(x$X1)
  p <- ncol(x$X)
  
  data <-cbind(x$X,x$y)
  if (x$nboot == 0){
    res <- estim.qosa(data, 1:n, type = x$type, alpha = x$alpha, n1 = n1, nsample = x$nsample)
    x$S <- data.frame(res)
    colnames(x$S) <- "original"
    rownames(x$S) <- colnames(x$X)
  }else{
    S.boot <- boot(data, estim.qosa, R = x$nboot, type = x$type, alpha = x$alpha, n1 = n1, nsample = x$nsample)
    x$S <- bootstats(S.boot, x$conf, "basic")
    rownames(x$S) <- colnames(x$X1)
  }
  
  assign(id, x, parent.frame())
}

print.qosa<- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (! is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    if (! is.null(x$S)) {
      if (x$type=="mean"){
        cat("\n\n\nFirst order Sobol indices: \n")
        print(x$S)
      }else{
        cat("\n\n\nFirst order quantile-oriented indices: \n")
        print(x$S)
      }
    }
    else{
      cat("(empty)\n")
    }
  }
}

plot.qosa <- function(x, ylim = c(0, 1), ...) {
  
  if (! is.null(x$y)) {
    if (x$type=="mean"){
      legend = "First-order Sobol indices"
    }else{
      legend = "First-order quantile-oriented indices"
    }
    nodeplot(x$S, ylim = ylim)
    legend(x = "topright", legend = legend)
  }
}

ggplot.qosa <- function(data, mapping = aes(), ylim = c(0,1), ..., environment = parent.frame()) {
  x <- data
  
  if (! is.null(x$y)) {
    if (x$type=="mean"){
      xname = "First-order Sobol indices"
    }else{
      xname = "First-order quantile-oriented indices"
    }
    nodeggplot(list(x$S), xname = xname, ylim = ylim)
  }
}