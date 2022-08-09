# library(boot)
#library(numbers)

# LAST MODIFIED:
# 16/12/2016

# AUTHORS:
# Gilquin Laurent

# SUMMARY:
# Extension of the replication procedure introduced by Tissot & Prieur (2015)
# to recursively estimate either first-order or closed second-order indices.

sobolrec <- function(model=NULL, factors, layers, order, precision, method=NULL, tail=TRUE, ...) {
  
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
  
  if (order==1 | order==2) {
    t <- order
  }
  else {
    stop("invalid argument 'order', waiting for value 1 or 2.")
  }
  
  if (order==2) {
    if (!(method=="ar" | method=="al")) {
      stop("invalid argument 'method', waiting for character 'al' or 'ar' ")
    }
  }  
  
  if(length(precision)!=2){
    stop("invalid argument 'precision', waiting for a vector of size 2")
  }
  
  if(precision[2]<2){
    stop("invalid argument 'precision', the second component must be greater than 2")
  }
  
  if (!is.logical(tail)) {
    stop("invalid argument 'tail', waiting for a boolean.")
  }
  
  # Conditions checking
  if (t==2) {
    if(layers>=(d-1)^2){
      if (sqrt(layers)%%1==0) {
        if (numbers::isPrime(sqrt(layers))){
          q <- sqrt(layers)
        }
        else if (length(unique(primeFactors(sqrt(layers))))==1){
          q <- sqrt(layers)
        } else {          
          if (tail){
            q <- numbers::previousPrime(sqrt(layers))
          } else {
            q <- numbers::nextPrime(sqrt(layers))
          }
          warning("The value entered for layers is not the square of a prime number. It has been replaced by: ",paste(q^2))
        }
      }
      if (sqrt(layers)%%1!=0) {
        if (tail){
          q <- numbers::previousPrime(sqrt(layers))
        } else {
          q <- numbers::nextPrime(sqrt(layers))
        }
        warning("The value entered for layers is not the square of a prime number. It has been replaced by: ",paste(q^2))
      }
    }  
    if(layers<d^2){
      q <- numbers::nextPrime(d)
      layers <- q^2
      warning("The value entered for layers is not satisfying the constraint. It has been replaced by: ",paste(q^2))
    }
  }
  
  # Construction of the nested Design of Experiments
  # case of strength 1   
  if (t==1) {
    n <- prod(layers)
    # main structure allocation
    doe <- matrix(NA,nrow=n,ncol=2*d)
    # matrix of permutations
    perm <- matrix(NA,nrow=n,ncol=2*d)
    perm[,1:d] <- nested_lhs_cplus(d,layers) 
    # matrix of random numbers
    mat_rand <- matrix(runif(d*n),nrow=n,ncol=d)
    # construction of the two nested designs
    layers <- c(0,cumprod(layers))
    for (i in 1:d) {
      for(j in 1:(length(layers)-1)){
        perm[(layers[j]+1):layers[j+1],d+i] <- sample(perm[(layers[j]+1):layers[j+1],i])
      }
      p1 <- as.numeric(perm[,c(i,d+i)])
      doe[,c(i,d+i)] <- (p1-mat_rand[p1,i])/n
    }
    conv <- matrix(1,nrow=d,ncol=length(layers)+1)
    p2 <- NULL
  }
  
  # case of strength 2
  if (t==2) {
    
    if(method=="al"){
      # number layers
      coeff <- min(floor(10^6/q^2),q^(d-2)) # NOTE: coeff is set to goes up to 10^6 points for each replicated designs
      # random choice of the blocks
      vect <- matrix(sample(0:(q-1),10^5*(d-2),replace=TRUE),ncol=d-2)
      vect <- t(unique(vect))
      if(d==3){
        vect <- matrix(vect[1:coeff],ncol=coeff)
      } else {
        vect <- vect[,1:coeff]
      }
      n <- q^2*coeff
      #construction of the nested orthogonal arrays
      OAa <- addelman_const(d,q,choice="U")-1
      OAb <- addelman_const(d,q,choice="W")-1
      gf_add <- field_tables(q)$add
      vect <- rbind(matrix(rep(0,ncol(vect)*2),nrow=2),vect)
      doe0a <- apply(outer(OAa,rep(0,ncol(vect)),'+'),2,I)
      doe0b <- apply(outer(OAb,rep(0,ncol(vect)),'+'),2,I)
      pairsa <- cbind(c(doe0a),rep(c(t(vect)),each=q^2))
      pairsb <- cbind(c(doe0b),rep(c(t(vect)),each=q^2))
      doe0a <- matrix(gf_add[pairsa+1],ncol=d)+1
      doe0b <- matrix(gf_add[pairsb+1],ncol=d)+1
      rm(vect,OAa,OAb,gf_add,pairsa,pairsb)
    }
    if (method=="ar"){
      #number layers
      coeff <- min(500,q^(d-2)*0.75) # NOTE: the accept-reject loop is times consuming, coeff is set to 500 max
      doe <- matrix(0,nrow=1,ncol=d)
      OAa <- addelman_const(d,q,choice="U")
      blok <- matrix(0,nrow=q^2,ncol=d)
      it <- 0
      count <- 0
      # recursive construction of the block
      while (it<(coeff) & count < 20){ # NOTE: count<20 is a time limit to break the loop if no more blocks are to be found in 20 consecutive trials
        perm <- replicate(d,sample(q))
        ind <- sample(d)
        for (i in 1:d){
          blok[,i] <- perm[OAa[,ind[i]],i]
        }
        bool <- !Compar_array(t(doe),t(blok))
        if (bool){
          doe <- rbind(doe,blok)
          it <- it + 1
          count <- 0
        } else {
          count <- count+1
        }
      }
      coeff <- it
      n <- q^2*coeff
      doe0a <- doe[-1,]
      doe0b <- doe0a
      rm(it,blok,doe,perm,ind,count)
    }
    # main structure allocation
    doe <- matrix(NA,nrow=n,ncol=2*d)
    # matrix of permutations
    perm <- replicate(2*d,sample(q))
    # matrix of random numbers
    mat_rand <- lhs(n,d)
    # construction of the two nested designs
    for (i in 1:d) {
      p1 <- perm[doe0a[,i],i]
      p2 <- perm[doe0b[,i],d+i]
      doe[,c(i,d+i)] <- (c(p1,p2)-mat_rand[c(p1+rep(q^2*seq(0,coeff-1),each=q^2),p2+rep(q^2*seq(0,coeff-1),each=q^2)),i])/q
    }
    layers <- c(0,outer(q^2,1:coeff,'*'))
    conv <- matrix(1,nrow=d*(d-1)/2,ncol=length(layers)+1)
  }
  
  #Stocking of the two replicated designs
  X <- rbind(doe[,1:d],doe[,(d+1):(2*d)])
  colnames(X) <- X.labels
  rm(doe,p1,p2,perm)
  
  # object of class "sobolrec"
  x <- list(model=model, factors=factors, X=X, layers=layers, order=t, precision=precision, method=method,
            conv=conv, call=match.call())
  class(x) <- "sobolrec"
  
  # computing the response if the model is given
  if (!is.null(x$model)) {
    k <- 1
    stop_crit <- FALSE
    while( !(stop_crit) & (k<length(layers))){
      ask(x, index=k, ...)
      tell(x, index=k, ...)
      stop_crit <- x$stop_crit
      k <- k+1
    }
  }
  return(x)
}

# ------------------------------------------------------------------
# Ask method to recover response and permutation matrix

ask.sobolrec <- function(x, index, ...){
  
  id <- deparse(substitute(x))
  d <- x$factors
  layers <- x$layers
  t <- x$order
  k <- index
  n <- nrow(x$X)/2
  if(t==1){
    loop_index <- matrix(1:d,ncol=1)
  } else {
    loop_index <- t(combn(d,2))
  }
  RP <- matrix(NA,nrow=layers[k+1]-layers[k],ncol=nrow(loop_index))
  block <- cbind(x$X[(layers[k]+1):layers[k+1],],x$X[n+(layers[k]+1):layers[k+1],])
  if(!is.null(x$model)){
    sobolrec.response(x, block_id=seq(layers[k]+1,layers[k+1]), ...)
  } else {
    x$block <- rbind(block[,1:d],block[,(d+1):(2*d)])
  }
  for(ind in 1:nrow(loop_index)){
    p1 <- do.call(base::order,as.data.frame(block[,loop_index[ind,]]))
    p2 <- do.call(base::order,as.data.frame(block[,d+loop_index[ind,]]))
    RP[,ind] <- p2[order(p1)]
  }
  x$RP <- RP
  assign(id, x, parent.frame())
}


# --------------------------------------------------------------------
# Estim method to estimate Sobol' indices 

estim.sobolrec <- function(data, i = 1 : nrow(data), RP, rec_val, ...){
  
  # local variables
  n <- nrow(data)
  Y_i <- data[i,1]
  nd <- ncol(RP)
  Y <- matrix(NA,ncol=nd,nrow=n)
  for(j in 1:nd){
    Y[,j] <- data[RP[i,j],2]
  }
  
  #Sobol' indices calculation through recursive formula
  N <- rec_val[length(rec_val)]
  phi <- (rec_val[1]+sum(Y_i))/N
  psi <- (rec_val[2:(nd+1)]+Y_i%*%Y)/N
  dzeta <- (rec_val[(nd+2):(2*nd+1)]+colSums(Y))/N
  V <- (rec_val[2*nd+2]+sum(Y_i^2)-N*phi^2)/N
  rec_val <- c(phi,psi,dzeta,V+phi^2)
  return(c(V,rec_val))
}


# --------------------------------------------------------------------
# Tell method to estimate Sobol' indices and compute bootstrap confidence intervals
 
tell.sobolrec <- function(x, y = NULL, index, ...){
  
  id <- deparse(substitute(x))
  if (is.null(x$S)) {
    x$S <-list()
  } 
  if (!is.null(y)) {
    x$y_out[[index]] <- y
    x$y <- y
  } 
  else if (is.null(x$y)) {
    stop("y not found")
  }
  
  # Sobol' indices estimation and confidence intervals
  RP <- x$RP
  d <- x$factors
  t <- x$order
  nd <- ncol(RP)
  data <- matrix(x$y,ncol=2)
  N <- x$layers[index+1]
  # retrieving values of recursive terms
  if(is.null(x$rec_val)){
    rec_val <- c(rep(0,2*nd+2),N)
  } else {
    rec_val <- c(x$rec_val,N)
  }
  # call to estim.sobolrec
  V <- data.frame(original = estim.sobolrec(data=data, RP=RP, rec_val=rec_val))
  S <- data.frame(original=(V[3:(nd+2),1]-V[2,1]*V[(nd+3):(2*nd+2),1]) / V[1,1])
  x$rec_val <- N*V[2:nrow(V),1]
  x$V[[index]] <- V[1,1] 
  
  # storage and stopping criterion test
  x$S[[index]] <- S
  l_0 <- x$precision[2]
  if(index > (l_0-1)){
    x$conv[,index+1] <- abs(((x$S[[index]])[,1]-(x$S[[index-1]])[,1]))
    x$stop_crit <- all(rowSums(x$conv[,(index-l_0+1):index]<x$precision[1])==l_0)
  } else {
    x$stop_crit <- FALSE
  }
  x$steps <- index
  x$N <- x$layers[index+1]
  
  # output purpose
  if(t==1){
    loop_index <- matrix(1:d,ncol=1)
  } else {
    loop_index <- t(combn(d,2))
  }
  rownames <- paste("X",apply(loop_index,1,paste,collapse=""),sep= "")
  rownames(x$S[[index]]) <- rownames
  
  assign(id, x, parent.frame())
}


# --------------------------------------------------------------------
# Print method to copy results: model variance, percentage of missing values and Sobol' estimates.

print.sobolrec <- function(x, ...) {  
  
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  cat("\nModel runs:", 2*x$N, "\n")
  if (!is.null(x$y)) {
    cat("\nModel variance:\n")
    print(unlist(x$V[x$steps]))
    if(x$order==1){
      cat("\nFirst-order indices:\n")
      print(data.frame(x$S[x$steps]))
    }
    if(x$order==2){
      cat("\nClosed second-order indices:\n")
      print(data.frame(x$S[x$steps]))
    }
  }
}


# --------------------------------------------------------------------
# Plot method to draw Sobol' estimates

plot.sobolrec <- function(x, ylim = c(0, 1), ...) {
  
  if (!is.null(x$y)) {
    if (x$order==1){
      nodeplot(x$S[[x$stopp]], ylim = ylim, ...)
      legend(x = "topright", legend = c("First-order indices"))
    }
    if (x$order==2){
      nodeplot(x$S2[[x$stopp]], ylim = ylim, ...)
      legend(x = "topright", legend = c("Closed second-order indices"))    
    }
  }
}


# response modified for sobolrec
sobolrec.response <- function(x, block_id = NULL, ...) {
  id <- deparse(substitute(x))
  
#  if (class(x$model) == "function") {
  if (inherits(x$model, "function")){
    if (!is.null(block_id)) {
      n <- nrow(x$X)
      y <- numeric(n)
      y <- x$model(x$X[c(block_id,n/2+block_id),], ...)
    } else {
      y <- x$model(x$X, ...)
    }
  } else if (TRUE %in% (paste("predict.", class(x$model), sep="") %in% methods(predict))) {
    y <- predict(x$model, x$X, ...)
  } else {
    stop("The model isn't a function or does not have a predict method")
  }
  
#  if (class(y) != "numeric") {
  if (!inherits(y, "numeric")){
    y <- as.numeric(y)
    warning("Conversion of the response to numeric")
  }
  
  x$y <- y
  assign(id, x, parent.frame())
}

# ----------------------------------------------------------------
# Function generating a LHS

lhs <- function (n, dimension){
  
  shift <- matrix(runif(n * dimension), nrow = n, ncol = dimension)
  perm <- replicate(dimension,sample(n))
  shift <- (perm-shift)/n
  return(shift)
}

