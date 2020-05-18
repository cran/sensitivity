# library(boot)
#library(numbers)

# LAST MODIFIED:
# 15/01/2018

# AUTHORS:
# Gilquin Laurent

# SUMMARY:
# Extension of the replication procedure introduced by Tissot & Prieur (2015)
# to estimate both first-order and second-order indices at a cost of 2*N.

sobolrep <- function(model=NULL, factors, N, tail=TRUE, conf=0.95, nboot=0, nbrep=1, total=FALSE, ...) {
  
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
    stop("invalid argument 'factors', expecting a positive integer or a character string vector (names).")
  }
  
  if (N%%1!=0 | N<=0) {
    stop("invalid argument 'N', expecting a positive integer.")
  }
  
  if (conf < 0 | conf > 1) {
    stop("invalid argument 'conf', expecting a value in ]0,1[.")
  }
  
  if (!is.logical(tail)) {
    stop("invalid argument 'tail', expecting a boolean.")
  }
  
  if(nboot%%1!=0 | nboot<0){
    stop("invalid argument 'nboot', expecting a positive integer or zero.")
  }
  
  if(!is.logical(total)){
    stop("invalid argument 'total', expecting a boolean.")
  }
  
  # Conditions checking:
  if(N>=(d-1)^2){
    if (sqrt(N)%%1==0) {
      if (numbers::isPrime(sqrt(N))){
        q <- sqrt(N)
      }
      else if (length(unique(primeFactors(sqrt(N))))==1){
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
    if (sqrt(N)%%1!=0) {
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
  if(N<d^2){
    q <- numbers::nextPrime(d)
    N <- q^2
    warning("The value entered for N is not satisfying the constraint N >= (d-1)^2. It has been replaced by: ",paste(q^2))
  }
  
  # Main structures allocation
  doe <- matrix(NA,nrow=N,ncol=2*d)
  # matrix of permutations
  perm <- replicate(2*d,sample(q))
  # orthogonal array
  doe0a <- addelman_const(d,q,choice="U")
  doe0b <- addelman_const(d,q,choice="W")
  # matrix of random numbers
  mat_rand <- matrix(runif(d*q),nrow=q,ncol=d)
  # construction of the two replicated designs
  for (i in 1:d) {
    p1 <- perm[doe0a[,i],i]
    p2 <- perm[doe0b[,i],d+i]
    doe[,c(i,d+i)] <- (c(p1,p2)-mat_rand[c(p1,p2),i])/q
  }
  
  # Stocking the ordering matrix
  loop_index <- t(combn(d,2))
  RP <- matrix(NA,nrow=N,ncol=d+nrow(loop_index))
  for(ind in 1:d){
    p1 <- order(doe[,ind])
    p2 <- order(doe[,d+ind])
    RP[,ind] <- p2[order(p1)]
  }
  for(ind in 1:nrow(loop_index)){
    p1 <- do.call(base::order,as.data.frame(doe[,loop_index[ind,]]))
    p2 <- do.call(base::order,as.data.frame(doe[,d+loop_index[ind,]]))
    RP[,d+ind] <- p2[order(p1)]
  }
  
  # Deleting unused objects
  rm(perm,doe0a,doe0b,mat_rand,p1,p2,loop_index)
  
  #Stocking of the two replicated designs
  X <- rbind(doe[,1:d],doe[,(d+1):(2*d)])
  
  # Dealing with total case
  if(total){
    
    # LHS sampling for the d additional designs
    perm <- replicate(d,sample(N))
    doetot <- matrix(runif(N*d),nrow=N,ncol=d)
    for (i in 1:d){
      p1 <- perm[,i]
      doetot[,i] <- (p1-doetot[p1,i])/N
    }
    
    # creation of the d addtional designs
    index <- seq(1,d) + seq(0,d-1)*d
    X0 <- matrix(c(rep(doe[,(d+1):(2*d)],d)),ncol=d^2)
    X0[,index] <- doetot
    Xtot <- matrix(X0[,c(outer(seq(1,d^2,by=d),seq(0,d-1),'+'))],ncol=d)
    X <- rbind(X, Xtot)
    rm(perm,doetot,p1,X0,index)
  }
  
  # object of class "sobolrep"
  x <- list(model=model, factors=factors, X=X, RP=RP, N=N, order=t,
            conf=conf, tail=tail, nboot=nboot, nbrep=nbrep, total=total, call=match.call())
  class(x) <- "sobolrep"
  
  # computing the response if the model is given
  if (!is.null(x$model)) { 
    response(x, ...) 
    tell(x, ...) 
  }
  return(x)
}

# --------------------------------------------------------------------
# Get new repetition

repeat.sobolrep <- function(Ord, RP, mat_rep = NULL, pow_init=NULL){
  
  new_RP <- matrix(0,ncol=ncol(RP),nrow=nrow(RP))
  d <- ncol(RP)
  q <- sqrt(nrow(RP))
  for(ind in 1:d){
    p1 <- RP[,ind]
    col_Ord <- matrix(Ord[,ind],ncol=q)
    new_col <- c(apply(col_Ord,2,sample))
    map <- setNames(new_col,c(col_Ord))
    new_RP[,ind] <- unname(map[as.character(p1)])
  }
  return(new_RP)
}

# --------------------------------------------------------------------
# Estim method to estimate first-order and second-order Sobol' indices 

#-------------------------
#  NEW CODE : Rcpp Version
#-------------------------

estim.sobolrep <- function(data, i = 1 : nrow(data), RP, d, I, ...){
  
  res <- cpp_get_indices(data, RP, I, i, d)
  out <- c(unlist(res))
  
  return(out)
}
  
  #-------------------------
  #  OLD CODE : R Version
  #-------------------------
  
# estim.sobolrep <- function(data, i = 1 : nrow(data), RP, d, I, ...){
#   
#   # local variables
#   N <- nrow(data)
#   Y_i <- data[i,1]
# 
#   #Sobol' indices calculation
#   
#   # first-order
#   Y <- matrix(NA,ncol=d,nrow=N)
#   for(j in 1:d){
#     Y[,j] <- data[RP[i,j],2]
#   }
# 
#   Mean <- colMeans(Y)
#   MeanY <- (mean(Y_i)+Mean)/2
# 
#   S <- rep(NA,d)
#   a <- rep(0,d)
#   b <- rep(0,d)
#   c <- rep(0,d)
#   out <- rep(0,d)
#   ind <- rep(1,d)
#   ind2 <- 2:(d+1)
#   S <- .C("LG_estimator",as.double(c(Y_i,Y)),as.double(MeanY),as.integer(d),as.integer(N),as.integer(ind),as.integer(ind2),as.double(a),as.double(b),as.double(c),as.double(out))[[10]]
# 
#   # second-order
#   Y <- matrix(NA,nrow=N,ncol=nrow(I))
#   for(j in 1:nrow(I)){
#     Y[,j] <- data[RP[i,d+j],2]
#   }
#   Mean <- colMeans(Y)
#   MeanY <- (mean(Y_i)+Mean)/2
# 
#   a <- rep(0,nrow(I))
#   b <- rep(0,nrow(I))
#   c <- rep(0,nrow(I))
#   out <- rep(0,nrow(I))
#   ind <- rep(1,nrow(I))
#   ind2 <- 2:(nrow(I)+1)
#   S2 <- .C("LG_estimator",as.double(c(Y_i,Y)),as.double(MeanY),as.integer(nrow(I)),as.integer(N),as.integer(ind),as.integer(ind2),as.double(a),as.double(b),as.double(c),as.double(out))[[10]]
#   S2 <- S2-S[I[,1]]-S[I[,2]]
#   return(c(S,S2))
# }

# --------------------------------------------------------------------
# Estim method to estimate total-effect Sobol' indices 

estimtotal.sobolrep <- function(data, i = 1 : nrow(data), ...){
  
  
  #total effect Sobol' indices
  ST <- cpp_get_total_indices(data, i)
  
  return(ST)
}


# --------------------------------------------------------------------
# Tell method to estimate Sobol' indices and compute bootstrap confidence intervals


tell.sobolrep <- function(x, y = NULL, ...){
  
  id <- deparse(substitute(x))
  if (! is.null(y)) {
    x$y <- y
  } 
  else if (is.null(x$y)) {
    stop("y not found")
  }
  
  # Arguments
  X2 <- x$X[(x$N+1):(2*x$N),]
  Ord <- apply(X2,2,base::order)
  RP <- x$RP
  d <- x$factors
  nbrep <- x$nbrep
  total <- x$total
  if(total){
    fulldata <- matrix(x$y,ncol=d+2)
    data <- fulldata[,1:2]
    datatot <- fulldata[,2:(d+2)]
  } else{
    data <- matrix(x$y,ncol=2)
  }
  
  # Repetitions
  S0 <- 0
  Sb <- 0
  S2b <- 0
  S20 <- 0
  timm <- as.numeric(Sys.time())
  seed_val <- ((timm-floor(timm))) * 1e8
  I <- t(combn(d,2))
  
  for(rep in 1:nbrep){
    new_RP <- cbind(repeat.sobolrep(Ord,RP[,1:d]),RP[,(d+1):ncol(RP)])
    # Sobol' indices estimation and confidence intervals
    if (x$nboot == 0){
      indices <- data.frame(original = estim.sobolrep(data=data, RP=new_RP, d=d, I=I))
      Sb <- data.frame("original"=indices[1:d,]) + Sb
      S2b <- data.frame("original"=indices[(d+1):nrow(indices),]) + S2b
    } else{
      set.seed(seed_val)
      S.boot <- boot(data=data, estim.sobolrep, R = x$nboot, RP=new_RP, d=d, I=I)
      Sb <- as.matrix(S.boot$t[,1:d]) + Sb
      S0 <- S.boot$t0[1:d] + S0
      S2b <- as.matrix(S.boot$t[,(d+1):ncol(RP)]) + S2b
      S20 <- S.boot$t0[(d+1):ncol(RP)] + S20
    }
    set.seed(NULL)
  }
  if (x$nboot == 0){
    x$S <- Sb/nbrep
    x$S2 <- S2b/nbrep
  } else{
    S.boot$t <- cbind(Sb, S2b)/x$nbrep
    S.boot$t0 <- c(S0, S20)/x$nbrep
    indices <- bootstats(S.boot, x$conf, "bias corrected")
    x$S <- indices[1:d,]
    x$S2 <- indices[(d+1):nrow(indices),]
  }
  if(total){
    if (x$nboot == 0){
      x$T <- data.frame("original"=estimtotal.sobolrep(datatot))
    } else{
      ST.boot <- boot(datatot, estimtotal.sobolrep, R = x$nboot)
      x$T <- bootstats(ST.boot, x$conf, "bias corrected")
    }
  }
  
  # output
  x$V <- var(data[,1])
  rownames <- paste("X",1:d,sep="")
  rownames2 <- paste("X",apply(t(combn(d,2)),1,paste,collapse=""),sep= "")
  rownames(x$S) <- rownames
  rownames(x$S2) <- rownames2
  if(total){
    rownames(x$T) <- rownames
  }
  
  assign(id, x, parent.frame())
}


# --------------------------------------------------------------------
# Print method to copy results: model variance, percentage of missing values and Sobol' estimates.

print.sobolrep <- function(x, ...) {  
  
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if(x$total){
    cat("\nModel runs:", x$N*(x$factors+2), "\n")
  } else{
    cat("\nModel runs:", 2*x$N, "\n")
  }
  if (!is.null(x$y)) {
    cat("\nModel variance:\n")
    print(x$V)
    cat("\nFirst-order indices:\n")
    print(x$S)
    cat("\nClosed second-order indices:\n")
    print(x$S2)
    if(x$total){
      cat("\nTotal-effect indices:\n")
      print(x$T)
    }
  }
}


# --------------------------------------------------------------------
# Plot method to draw Sobol' estimates

plot.sobolrep <- function(x, ylim = c(0, 1), choice, ...) {
  
  if (!is.null(x$y)) {
    if (choice==1){
      nodeplot(x$S, ylim = ylim, ...)
      legend(x = "topright", legend = c("First-order indices"))
    }
    if (choice==2){
      nodeplot(x$S2, ylim = ylim, ...)
      legend(x = "topright", legend = c("Second-order indices"))    
    }
    if (choice==3 & x$total){
      nodeplot(x$T, ylim = ylim, ...)
      legend(x = "topright", legend = c("Total-effect indices"))    
    }
  }
}
