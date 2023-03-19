# library(boot)
# library(numbers)

# Sub-files:
# * sobolroa_subroutines.R

# LAST MODIFIED: 
# December 9, 2016

# AUTHOR:
# Gilquin Laurent

# SUMMARY:
# Extension of the replication procedure (Tissot and Prieur 2015) to deal with sets of constraints (Gilquin et al., 2015). 
# This procedure estimates either all first-order or all second-order Sobol indices given two
# replicated designs.

sobolroauc=function(model=NULL, factors, constraints=NULL, N, p=1, order, tail=TRUE, conf=0.95, nboot=0, ...) {
  
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
  
  if(!is.list(constraints)) {
    stop("invalid argument 'constraints', waiting for a list")
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
  
  # factors resizing 
  const <- length(unlist(constraints)) # number of constrained variables
  const_set <- length(constraints) # number of set of constrained variables
  newd <- d-const+const_set # number of set variables
  if(const_set!=0){
    ind <- (1:d)[-unlist(constraints)]
    ind_length <- length(ind)
  } else {
    ind <- 1:d
    ind_length <- length(ind)
  }
  
  
  # Conditions checking for case order==2:
  if (t==2) {
    
    if(N>=(newd-1)^2){
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
    
    if(N<(newd-1)^2){
      
      if (requireNamespace("numbers", quietly = TRUE)){
        if(numbers::isPrime(newd-1)){
          q <- newd-1
          N <- q^2
        } else {
          q <- numbers::nextPrime(newd-1)
          N <- q^2
        }
      }
      warning("The value entered for N is not satisfying the constraint N >= (newd-1)^2. It has been replaced by: ",paste(q^2))
    }
  }

  # Main structure allocation
  doe <- matrix(NA,nrow=N,ncol=2*d)
  if(t==1){
    # matrix of permutations
    perm <- replicate(2*newd,sample(N))
    # for standardisation 
    doe0 <- replicate(newd,seq(1,N))
    # matrix of random numbers
    mat_rand <- matrix(runif(d*N),nrow=N,ncol=d)
  } else {
    # matrix of permutations
    perm <- replicate(2*newd,sample(q))
    # orthogonal array
    doe0 <- apply(outer(sobolroa.Hadam(q)[,1:newd],0:(q-1),'+')%%q+1,2,I)
    if (q==newd-1){
      doe0 <- cbind(doe0,rep(seq(1,q),q))
    }
    # matrix of random numbers
    mat_rand <- matrix(runif(d*q),nrow=q,ncol=d)
  }
  
  # construction of the two replicated designs
  discret <- ifelse(t==1,N,q)
  # sampling of the non-constrained factors
  if (ind_length!=0){
    for (i in 1:ind_length){
      ind_col <- ind[i]
      p1 <- perm[doe0[,i],c(i,d+i)]
      doe[,c(ind_col,d+ind_col)] <- (p1-mat_rand[p1,ind_col])/discret
    }
  }
  # sampling of the constrained factors
  if (const_set!=0){
    for (i in 1:const_set){
      numb_col <- d-const+i
      index <- constraints[[i]]
      index_length <- length(index)
      # call to rowsort (external C function)
      C_x <- c(t(mat_rand[,index]))
      C_inter <- rep(0,index_length)
      C_out <- rep(0,length(C_x))
      C_res <- .C("LG_rowsort",as.integer(discret),as.double(C_x),as.integer(index_length),
                  as.double(C_inter),as.double(C_out))[[5]]
      #Space-filling simplex sampling
      U <- matrix(C_res,ncol=index_length,byrow=TRUE)
      U <- c(t(matrix(c(rep(0,discret),U,rep(1,discret)),ncol=index_length+2)))
      S <- U[-1]-U[-length(U)]
      S <- matrix(S[S>0],ncol=index_length+1,byrow=TRUE)
      #Adjust the subdivision
      test <- length(unlist(strsplit(paste(discret^(1/index_length)),"")))<5
      if (!test){
        numb_simplex <- floor(discret^(1/index_length))+1
        simplex <- sobolroauc.simplex_create(numb_simplex+1,index_length)
        drop <- sample(1:numb_simplex^index_length,numb_simplex^index_length-discret)
        P <- matrix(c(t(S))*simplex$p[t(simplex$s[-drop,]),],nrow=index_length+1)
      } else {
        numb_simplex <- floor(discret^(1/index_length)-10^(-5))+1
        simplex <- sobolroauc.simplex_create(numb_simplex+1,index_length)
        P <- matrix(c(t(S))*simplex$p[t(simplex$s),],nrow=index_length+1)
      }
      #Fill the design
      P <- matrix(colSums(P),ncol=index_length)
      doe[,index] <- P[perm[doe0[,numb_col],numb_col],]
      doe[,d+index] <- P[perm[doe0[,numb_col],newd+numb_col],]
    }
    #Deleting unused objects
    rm(C_x,C_inter,C_out,U,S,P,numb_simplex,simplex,drop)
  }
    
  #Stocking the ordering matrix
  if(t==1){
    loop_index <- matrix(1:newd,ncol=1)
  } else {
    loop_index <- t(utils::combn(newd,2))
  }
  RP <- matrix(NA,nrow=N,ncol=nrow(loop_index))
  if(const_set!=0){
    col_ind <- c(ind,sapply(constraints, "[[", 1))
  } else {
    col_ind <- 1:d
  }
  for(ind in 1:nrow(loop_index)){
    p1 <- do.call(base::order,as.data.frame(doe[,col_ind[loop_index[ind,]]]))
    p2 <- do.call(base::order,as.data.frame(doe[,d+col_ind[loop_index[ind,]]]))
    RP[,ind] <- p2[order(p1)]
  }
  
  #Deleting unused objects
  rm(perm,mat_rand,p1,p2,loop_index,discret,col_ind)
  if(t==1){
    doe0 <- NULL
  }
  
  #Stocking of the two replicated designs
  X <- rbind(doe[,1:d],doe[,(d+1):(2*d)])
  colnames(X) <- X.labels
  
  # object of class "sobolroauc"
  x <- list(model=model, factors=factors, constraints=constraints, X=X, OA=doe0, N=N, p=p, order=t,
            conf=conf, tail=tail, nboot=nboot, RP=RP, call=match.call())
  class(x) <- "sobolroauc"
  
  # computing the response if the model is given
  if (!is.null(x$model)) { 
    response(x, ...) 
    tell(x, ...) 
  }
  return(x)
}


# --------------------------------------------------------------------
# Estim method to estimate Sobol' indices 


estim.sobolroauc=function(data, i = 1 : nrow(data), RP, t, p, ...){
  
  # local variables
  nd <- ncol(RP)
  newd <- ifelse(t==1,nd,Re(polyroot(c(-2*nd,-1,1))[1]))
  N <- nrow(RP)
  Ya <- matrix(data[i,seq(1,2*p,by=2)],ncol=p)
  Yb <- matrix(data[,seq(1,2*p,by=2)+1],ncol=p)
  
  
  #Missing values purpose
  na.rm <- TRUE
  Na <- apply(Ya,2,function(x){as.numeric(sum(!is.na(x)))})
  
  #Variance calculation
  V <- sum(diag(var(Ya, na.rm=TRUE)))
  
  # Loop index
  if(t==1){
    loop_index <- matrix(1:newd,ncol=1)
  } else {
    loop_index <- t(utils::combn(newd,2))
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


tell.sobolroauc=function(x, y = NULL, ...){
  
  id <- deparse(substitute(x))
  if (! is.null(y)) {
    x$y <- y
  } 
  else if (is.null(x$y)) {
    stop("y not found")
  }
  
  # Sobol' indices estimation and confidence intervals
  d <- x$factors
  constraints <- x$constraints
  t <- x$order
  p <- x$p
  RP <- x$RP
  nd <- ncol(RP)
  data <- matrix(x$y,nrow=x$N)
  
  if (x$nboot == 0){
    V <- data.frame(original = estim.sobolroauc(data=data, RP=RP, t=t, p=p))
    S <- V[2:(nd + 1), 1, drop=FALSE] / V[1,1]
  } else{
    V.boot <- boot(data=data, estim.sobolroauc, R = x$nboot, RP=RP, t=t, p=p)
    V <- bootstats(V.boot, x$conf, "basic")
    S.boot <- V.boot
    S.boot$t0 <- V.boot$t0[2:(nd + 1)] / V.boot$t0[1]
    S.boot$t <- V.boot$t[,2:(nd + 1)] / V.boot$t[,1]
    S <- bootstats(S.boot, x$conf, "basic")
  }
  
  # output
  x$V <- V
  x$S <- S
  
  const <- length(unlist(constraints)) 
  const_set <- length(constraints) 
  newd <- d-const+const_set 
  if(t==1){
    loop_index <- matrix(1:newd,ncol=1)
  } else {
    loop_index <- t(combn(newd,2))
  }
  if (const_set!=0){
    name_list <- unlist(lapply(c((1:d)[-unlist(constraints)],constraints),function(x){paste("{",paste(x,collapse=","),"}",sep="")}))
  } else {
    name_list <- paste("{",1:d,"}",sep="")
  }
  loop_index <- matrix(name_list[loop_index],ncol=t)
  rownames <- paste("S",apply(loop_index,1,paste,collapse=""),sep= "")
  rownames(x$S) <- rownames
  rownames(x$V) <- c("global",gsub("S","V",rownames))
  
  assign(id, x, parent.frame())
}

# --------------------------------------------------------------------
# Print method to copy results: model variance, percentage of missing values and Sobol' estimates.

print.sobolroauc <- function(x, ...) {  
  
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

plot.sobolroauc <- function(x, ylim = c(0, 1), ...) {
  
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

ggplot.sobolroauc <- function(data, mapping = aes(), ylim = c(0,1), ..., environment = parent.frame()) {
  x <- data
  
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
