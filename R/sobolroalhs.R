# library(boot)
#library(numbers)

# Sub-files:
# * sobolroa_subroutines.R

# LAST MODIFIED:
# 09/12/2016

# AUTHORS:
# Gilquin Laurent
# Viry Laurence

# SUMMARY:
# Implementation of the replication procedure introduced by Tissot and Prieur (2015)

sobolroalhs=function(model=NULL, factors, N, p=1, order, tail=TRUE, conf=0.95, nboot=0, ...) {
  
  #Initialisation
  
  if (is.character(factors)) {
    X.labels <- factors
    d <- length(X.labels)
  }
  else if (factors%%1==0 & factors>0) {
    d <- factors
    X.labels <- paste("X", 1:d, sep = "")
  }
  else {
    stop("invalid argument 'factors', waiting for a positive integer or a character string vector (names).")
  }
  
  if (N%%1!=0 | N<=0) {
    stop("invalid argument 'N', waiting for a positive integer.")
  }
  
  if (p%%1!=0 | p<=0) {
    stop("invalid argument 'p', waiting for a positive integer.")
  }
  
  if (order!=1 | order!=2) {
    t <- order
  }
  else {
    stop("invalid argument 'order', waiting for value 1 or 2.")
  }
  
  if (conf < 0 | conf > 1) {
    stop("invalid argument 'conf', waiting for a value in ]0,1[.")
  }
  
  if (!is.logical(tail)) {
    stop("invalid argument 'tail', waiting for a boolean.")
  }
  
  if(nboot%%1!=0 | nboot<0){
    stop("invalid argument 'nboot', waiting for a positive integer or zero.")
  }
  
  # Conditions checking for case order==2:
  if (t==2) {
    
    if(N>=(d-1)^2){
      q <- sqrt(N)
      
      if (sqrt(N)%%1==0) {
        
        if (requireNamespace("numbers", quietly = TRUE)){ 
          
          if (numbers::isPrime(sqrt(N))){
            q <- sqrt(N)
          } else {          
            if (tail){
              q <- numbers::previousPrime(sqrt(N))
              N <- q^2
            } else {
              q <- numbers::nextPrime(sqrt(N))
              N <- q^2
            }
            warning("The value entered for N is not the square of a prime number. It has been replaced by: ",paste(q^2))
          }
        }
      }
      
      if (sqrt(N)%%1!=0) {
        if (requireNamespace("numbers", quietly = TRUE)){
          if (tail){
            q <- numbers::previousPrime(sqrt(N))
            N <- q^2
          } else {
            q <- numbers::nextPrime(sqrt(N))
            N <- q^2
          }
        }  
        warning("The value entered for N is not the square of a prime number. It has been replaced by: ",paste(q^2))
      }
    }  
    
    if(N<(d-1)^2){
      
      if (requireNamespace("numbers", quietly = TRUE)){
        if(numbers::isPrime(d-1)){
          q <- d-1
          N <- q^2
        } else {
          q <- numbers::nextPrime(d-1)
          N <- q^2
        }
      }
      warning("The value entered for N is not satisfying the constraint N >= (d-1)^2. It has been replaced by: ",paste(q^2))
    }
  }
  
  # Main structures allocation
  doe <- matrix(NA,nrow=N,ncol=2*d)
  if(t==1){
    # matrix of permutations
    perm <- replicate(2*d,sample(N))
    # for standardisation 
    doe0 <- replicate(d,seq(1,N))
    # matrix of random numbers
    mat_rand <- matrix(runif(d*N),nrow=N,ncol=d)
  } else {
    # matrix of permutations
    perm <- replicate(2*d,sample(q))
    # orthogonal array
    doe0 <- apply(outer(sobolroa.Hadam(q)[,1:d],0:(q-1),'+')%%q+1,2,I)
    if (q==d-1){
      doe0 <- cbind(doe0,rep(seq(1,q),q))
    }
    # matrix of random numbers
    mat_rand <- matrix(runif(d*q),nrow=q,ncol=d)
  }
  
  # construction of the two replicated designs
  discret <- ifelse(t==1,N,q)
  for (i in 1:d) {
    p1 <- perm[doe0[,i],c(i,d+i)]
    doe[,c(i,d+i)] <- (p1-mat_rand[p1,i])/discret
  }
  
  # Stocking the ordering matrix
  if(t==1){
    loop_index <- matrix(1:d,ncol=1)
  } else {
    loop_index <- t(combn(d,2))
  }
  RP <- matrix(NA,nrow=N,ncol=nrow(loop_index))
  for(ind in 1:nrow(loop_index)){
    p1 <- do.call(base::order,as.data.frame(doe[,loop_index[ind,]]))
    p2 <- do.call(base::order,as.data.frame(doe[,d+loop_index[ind,]]))
    RP[,ind] <- p2[order(p1)]
  }
  
  # Deleting unused objects
  rm(perm,mat_rand,p1,p2,loop_index,discret)
  if(t==1){
    doe0 <- NULL
  }
  
  #Stocking of the two replicated designs
  X <- rbind(doe[,1:d],doe[,(d+1):(2*d)])
  colnames(X) <- X.labels
  
  # object of class "sobolroalhs"
  x <- list(model=model, factors=factors, X=X, OA=doe0, N=N, p=p, order=t,
            conf=conf, tail=tail, nboot=nboot, RP=RP, call=match.call())
  class(x) <- "sobolroalhs"
  
  # computing the response if the model is given
  if (!is.null(x$model)) { 
    response(x, ...) 
    tell(x, ...) 
  }
  return(x)
}


# --------------------------------------------------------------------
# Estim method to estimate Sobol' indices 


estim.sobolroalhs=function(data, i = 1 : nrow(data), RP, t, p, ...){
  
  # local variables
  nd <- ncol(RP)
  d <- ifelse(t==1,nd,Re(polyroot(c(-2*nd,-1,1))[1]))
  N <- nrow(RP)
  Ya <- matrix(data[i,seq(1,2*p,by=2)],ncol=p)
  Yb <- matrix(data[,seq(1,2*p,by=2)+1],ncol=p)

  
  # for missing values purpose
  na.rm <- TRUE
  Na <- apply(Ya,2,function(x){as.numeric(sum(!is.na(x)))})
  
  #Variance calculation
  V <- sum(diag(var(Ya, na.rm=TRUE)))
  
  # Loop index
  if(t==1){
    loop_index <- matrix(1:d,ncol=1)
  } else {
    loop_index <- t(combn(d,2))
  }
  
  #Sobol' indices calculation 
  VCE <- rep(NA,ncol=nrow(loop_index))
  for (ind in 1:nrow(loop_index)) {
      Y1 <- Ya
      Y2 <- matrix(Yb[RP[i,ind],],ncol=p)
      Nb <- apply(Y2,2,function(x){as.numeric(sum(!is.na(x)))})
      Nab <- apply(Y1+Y2,2,function(x){as.numeric(sum(!is.na(x)))})
      VCE[ind] <- sum(colSums(Y1*Y2,na.rm=na.rm)/Nab-colSums(Y1,na.rm=na.rm)*colSums(Y2,na.rm=na.rm)/(Na*Nb))
  }
  
  return(c(V,VCE))
}

# --------------------------------------------------------------------
# Tell method to estimate Sobol' indices and compute bootstrap confidence intervals


tell.sobolroalhs=function(x, y = NULL, ...){
  
  id <- deparse(substitute(x))
  if (! is.null(y)) {
    x$y <- y
  } 
  else if (is.null(x$y)) {
    stop("y not found")
  }
  
  # Sobol' indices estimation and confidence intervals
  d <- x$factors
  t <- x$order
  p <- x$p
  RP <- x$RP
  nd <- ncol(RP)
  data <- matrix(x$y,nrow=x$N)
  
  if (x$nboot == 0){
    V <- data.frame(original = estim.sobolroalhs(data=data, RP=RP, t=t, p=p))
    S <- V[2:(nd + 1), 1, drop=FALSE] / V[1,1]
  } else{
    V.boot <- boot(data=data, estim.sobolroalhs, R = x$nboot, RP=RP, t=t, p=p)
    V <- bootstats(V.boot, x$conf, "basic")
    S.boot <- V.boot
    S.boot$t0 <- V.boot$t0[2:(nd + 1)] / V.boot$t0[1]
    S.boot$t <- V.boot$t[,2:(nd + 1)] / V.boot$t[,1]
    S <- bootstats(S.boot, x$conf, "basic")
  }
  
  # output
  x$V <- V
  x$S <- S
  
  if(x$order==1){
    loop_index <- matrix(1:d,ncol=1)
    rownames <- paste("X",apply(loop_index,1,paste,collapse=""),sep= "")

  } else {
    loop_index <- t(combn(d,2))
    rownames <- paste("X",apply(loop_index,1,paste,collapse="X"),sep= "")
  }
  
  rownames(x$S) <- rownames
  rownames(x$V) <- c("global",gsub("S","V",rownames))
  
  assign(id, x, parent.frame())
}


# --------------------------------------------------------------------
# Print method to copy results: model variance, percentage of missing values and Sobol' estimates.

print.sobolroalhs <- function(x, ...) {  
  
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  cat("\nModel runs:", 2*x$N, "\n")
  if (!is.null(x$y)) {
    if (any(is.na(x$y))){
      tauxNa <- as.numeric(length(x$y[is.na(x$y)]))/length(x$y)*100
      cat("\nPercentage of missing values(NA):", tauxNa, "\n")
    }
    cat("\nModel variance:\n")
    print(x$V[1,1])
    if (x$order==1){
      cat("\nFirst-order indices:\n")
      print(x$S)
    }
    if (x$order==2){
      cat("\nClosed second-order indices:\n")
      print(x$S)
    }
  }
}

# --------------------------------------------------------------------
# Plot method to draw Sobol' estimates

plot.sobolroalhs <- function(x, ylim = c(0, 1), ...) {
  
  if (!is.null(x$y)) {
    nodeplot(x$S, ylim = ylim, ...)
    if (x$order==1){
      legend(x = "topright", legend = c("First-order indices"))
    }
    else{
      legend(x = "topright", legend = c("Closed second-order indices"))    
    }
  }
}

ggplot.sobolroalhs <- function(x, ylim = c(0, 1), ...) {
  
  if (!is.null(x$y)) {
    if (x$order==1){
      title <- "First-order indices"
    }
    else{
      title <- "Closed second-order indices"   
    }
    nodeggplot(listx = list(x$S), xname = title, ylim = ylim, title = title)
  }
}
