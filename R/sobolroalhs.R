# Implementation of the Sobol's method version introduced by Tissot & Prieur. 
# This method estimate all first or second order Sobol indices via two replicated latin hypercubes.
# Author: Laurent Gilquin, INRIA, France
#
# imputs:
# q : number of levels of orthogonal array
# t : strength of orthogonal array
# alpha : risk of the asymptotical confidence interval

# outputs:
# S : vector (resp. 2 dimensional matrix) containing the d (resp. d*(d-1)/2) estimations of 
#          first order (resp. second order) Sobol indices
# E : vector (resp. 2 dimensional matrix) containing the d (resp. d*(d-1)/2) estimations of 
#         the radius of the confidence interval for the first order (resp. second order) Sobol indices
# V : total variance estimation


#GOOGLE SCHOOLAR : Heydayat & Al. 1999 orthogonal array

sobolroalhs=function(model=NULL, factors, levels, order, choice="A", conf=0.95,...) {

  #Initialisation
  
  if (is.character(factors)) {
    X.labels <- factors
    d <- length(X.labels)
  }
  else if (is.numeric(factors)) {
    d <- factors
    X.labels <- paste("X", 1:d, sep = "")
  }
  else {
    stop("invalid argument 'factors', waiting for a scalar (number) or a character string vector (names)")
  }
    
  if (is.numeric(levels)) {
    q <- levels
  }
  else {
    stop("invalid argument 'levels', waiting for a scalar (number) ")
  }
  
  if (order==1 | order==2) {
    t <- order
  }
  else {
    stop("invalid argument 'order', waiting for value 1 or 2 ")
  }
  
  if (choice=='A' | choice=='B') {
    c <- choice
  }
  else {
    stop("invalid argument 'choice', waiting for value A or B") 
  }
  
  if (conf < 0 | conf > 1) {
    stop("invalid argument 'conf', waiting for a value in ]0,1[")
  }
  
  # Construction of the Design of Experiments
  # case of strength 1   
  
  if (t==1) {
    n <- q
    # construction of two random matrix of permutations
    permu1 <- (replicate(d,sample(q))) 
    permu2 <- (replicate(d,sample(q)))
    # construction of the two design
    doe1 <- doe2 <- matrix(nrow=q,ncol=d)
    # construction of the random matrix 
    shift <- matrix(runif(d*q),nrow=q,ncol=d) 
    # construction of the two rlh
    for (i in 1:d) {
      p1 <- permu1[,i]
      p2 <- permu2[,i]
      doe1[,i] <- (p1-shift[p1,i])/q
      doe2[,i] <- (p2-shift[p2,i])/q
    }
    doe0 <- NULL # for ouput purpose 
    doe <- rbind(doe1, doe2)
  }
    
  # case of strength 2
  if (t==2) {
    if (!isPrime(levels)) {
      q <- nextPrime(sqrt(levels))
      warning("the number of levels recquired was not a prime number, the number was replaced by : ",paste(q))
    }
    if(q<(d-1)^2){
      q = (d-1)^2 
      warning("the number of levels recquired was not satisfying the constrain, the number was replaced by : ",paste(q))
    }
    n <- q^2
    # initialisation of the OA of strength 2
    doe0 <- matrix(nrow=n,ncol=d)
    # construction of two random matrix of permutations
    permu1 <- (replicate(d,sample(q)))
    permu2 <- (replicate(d,sample(q)))
    #construction of the OA of strength 2 based on Bose's construction
    if (c == "A"){
      for (i in 1:d) {
        doe0[,i] <- matrix(1+(repmat(matrix(seq(i,q*i,i),nrow=1),q,1)+repmat(matrix(seq(1,q),ncol=1),1,q))%%q,nrow=1,ncol=n)
      }
    }
    if (c == "B") {
      H <- Hadam(q)
      s <- matrix(seq(0,(q-1)))
      doe0 <- t(matrix(apply(s,1,'+',H[(1:d),]),ncol=q^2)%%q)+1
    }
    doe1 <- doe2 <- matrix(nrow=n,ncol=d)
    shift <- matrix(runif(d*q),nrow=q,ncol=d)
    #construction of the two oalh
    for (i in 1:d) {
      p1 <- permu1[doe0[,i],i]
      p2 <- permu2[doe0[,i],i]
      doe1[,i] <- (p1-shift[p1,i])/q
      doe2[,i] <- (p2-shift[p2,i])/q
    }
    doe <- rbind(doe1, doe2)
  }

  #Stocking DoE 
  X <- matrix(doe,ncol=d)
  colnames(X) <- X.labels
    
  # object of class "sobolroalhs"
  x <- list(model=model, factors=factors, X=X, levels=q, size=n, order=t, choice=c,
  conf=conf, permu1=permu1, permu2=permu2, doe0=doe0, call=match.call())
  class(x) <- "sobolroalhs"

  # computing the response if the model is given
  if (!is.null(x$model)) {
    response(x, ...)
    tell(x)
  }

  return(x)
  }



estim.sobolroalhs=function(x,choice){
  
  Y <- x$y
  X <- x$X
  d <- x$factors
  q <- x$levels
  t <- x$order
  n <- x$size
  alpha <- 1-x$conf
  permu1 <- x$permu1
  permu2 <- x$permu2
  doe0 <- x$doe0
  
  Ya <- Y[1:n]
  Yb <- Y[(n+1):(2*n)]
  m1 <- mean(Ya)
  m2 <- mean(Yb)
  Y1 <- Y2 <- rep(NA,n)

  if (choice=="variance"){
   Vs <- 1/n*sum(Ya^2)-(1/n*sum(Ya))^2 # standard total variance estimator
   Veff <- 1/(2*n)*sum(Ya^2+Yb^2)-(1/(2*n)*sum(Ya+Yb))^2 # Monod & Al. total variance estimator
   L=c(Vs,Veff)
  }
  
  if (choice=="indices"){
   # case of strength 1
   if (t==1) {
     S <- Seff <- E <- Eff <- matrix(nrow=1,ncol=d)
     for (i in 1:d) {
       Y1[permu1[,i]] <- Ya
       Y2[permu2[,i]] <- Yb
       S[i] <- (1/n*sum(Y1*Y2)-1/n^2*sum(Y1)*sum(Y2))/x$Vs  #standard sobol indices estimator
       Seff[i] <- (1/n*sum(Y1*Y2)-(1/(2*n)*sum(Y1+Y2))^2)/x$Veff #Monod & Al. sobol indices estimator
       E[i] <- (std((Y1-m1)*(Y2-m2)-S[i]/2*((Y1-m1)^2+(Y2-m2)^2))/std(Y))^2/sqrt(n)*qnorm(1-alpha/2)
       Eff[i] <- (std((Y1-m1)*(Y2-m2)-Seff[i]/2*((Y1-m1)^2+(Y2-m2)^2))/std(Y))^2/sqrt(n)*qnorm(1-alpha/2)
     }
   }  
    
   # case of strength 2
   if (t==2) {
     S <- Seff <- E <- Eff <- rep(NA,d*(d-1)/2)
     i <- 1
     j <- 1
     for (l in 1:(d*(d-1)/2)) {
       i <- i+(j==d)
       j <- ((j==d)*i + (j<d)*j)+1
       p1i <- permu1[doe0[,i],i]
       p1j <- permu1[doe0[,j],j]
       p2i <- permu2[doe0[,i],i]
       p2j <- permu2[doe0[,j],j]
       Y1[p1i + p1j*q-q] <- Ya
       Y2[p2i + p2j*q-q] <- Yb
       S[l] <- (1/n*sum(Y1*Y2)-1/n^2*sum(Y1)*sum(Y2))/x$Vs  #standard sobol indices estimator
       Seff[l] <- (1/n*sum(Y1*Y2)-(1/(2*n)*sum(Y1+Y2))^2)/x$Veff #Monod & Al. sobol indices estimator
       E[l] <- (std((Y1-m1)*(Y2-m2)-S[l]/2*((Y1-m1)^2+(Y2-m2)^2))/std(Y))^2/sqrt(n)*qnorm(1-alpha/2)
       Eff[l] <- (std((Y1-m1)*(Y2-m2)-Seff[l]/2*((Y1-m1)^2+(Y2-m2)^2))/std(Y))^2/sqrt(n)*qnorm(1-alpha/2)
     }
   }
   L=c(S,Seff,S-E,Seff-Eff,S+E,Seff+Eff)
  }
  return(L)

}  
 
  


tell.sobolroalhs=function(x, y = NULL, ...){
  
  id <- deparse(substitute(x))
  if (! is.null(y)) {
    x$y <- y
  } 
  else if (is.null(x$y)) {
    stop("y not found")
  }
  
  V <- estim.sobolroalhs(x,choice="variance")
  x$V <- data.frame(original=V)
  x$Vs <- V[1]
  x$Veff <- V[2]
  
  lab <- c("original", "min. c.i.", "max. c.i.")
  M <- matrix(estim.sobolroalhs(x, choice="indices"),ncol=3,dimnames=list(NULL,lab))
  x$S <- as.data.frame(round(M,6))  
  
  rownames=c()
  d=x$factors
  if (x$order==1){
    for (i in 1:d) {
      rownames=c(rownames,paste("S",i,sep=""))
    }
  }
  if (x$order==2){
    i <- 1
    j <- 1
    for (l in 1:(d*(d-1)/2)) {
      i <- i+(j==d)
      j <- ((j==d)*i + (j<d)*j)+1
      rownames=c(rownames,paste("S",i,j,sep=""))
    }
  }
  rownames(x$S)=c(rownames,gsub("S","Seff",rownames))
  
  x[c(5,9,10,11,15,16)] <- NULL
  assign(id, x, parent.frame())
}



print.sobolroalhs <- function(x, ...) {
  
  rownames=c()
  d=x$factors
  
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  cat("\nModel runs:", length(x$y), "\n")
  
  if (!is.null(x$y)) {
    cat("\nModel variance:\n")
    rownames(x$V)=c("Vs","Veff")
    print(x$V)
    if (x$order==1){
      cat("\nFirst order indices:\n")
      print(x$S)
    }
    if (x$order==2){
      cat("\nGeneralized second order indices:\n")
      print(x$S)
    }
  }
}



plot.sobolroalhs <- function(x, ylim = c(0, 1), type="standard", ...) {
  if (!is.null(x$y)) {
    l <- nrow(x$S)/2
    if(type=="standard"){
      data <- x$S[1:l,]
    }
    if(type=="monod"){
      data <- x$S[(l+1):(2*l),]
    }
    if(type!="standard" & type!="monod"){
      stop("invalid argument 'type'")
    }
    nodeplot(data, ylim = ylim, ...)
    if (x$order==1){
      legend(x = "topright", legend = c("First order indices"))
    }
    else{
      legend(x = "topright", legend = c("Second order subset indices"))      
    }
  }
}

# COMPLEMENTARY FUNCTION


# empirical variance

std <-function(x){
  n <- length(x);
  return(sqrt((n-1)/n)*sd(x))
}

# repmat matlab function

repmat <- function(X,nr,nc){
  if (is.matrix(X)) {
    rx <- dim(X)[1]
    cx <- dim(X)[2]
    M <- matrix(t(matrix(X,rx,cx*nc)),rx*nr,cx*nc,byrow=T)
  }
  else {
    stop("invalid argument 'X', waiting for a matrix ")
  }
  return(M)
}

# Hadamard matrix

Hadam=function(q){
  S <- seq(0,(q-1))
  return(S%*%t(S)%%q)
}
