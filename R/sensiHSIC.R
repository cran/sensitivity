##################################################################
# Sensitivity indices based on kernel embeddings of distributions
#
# Sebastien Da Veiga 2014, Anouar Meynaoui & Amandine Marrel 2019
##################################################################

# Kernel functions
rbf_hsic <- function(x,param,d=NULL){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  d <- d/param
  return(exp(-0.5*d^2))
}

laplace_hsic <- function(x,param,d=NULL){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  d <- d/param
  return(exp(-d))
}

dcov_hsic <- function(x,param,d=NULL){
  nobs <- length(x)
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  return(0.5*(-d^param+matrix(abs(x),nobs,nobs)^param+t(matrix(abs(x),nobs,nobs))^param))
}

raquad_hsic <- function(x,param,d=NULL){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  d <- d/param
  return(1-d^2/(1+d^2))
}

invmultiquad_hsic <- function(x,param,d=NULL){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  d <- d/param
  return(1/sqrt(d^2+1))
}

linear_hsic <- function(x,...){
  return(x%*%t(x))
}

matern3_hsic <- function(x,param,d=NULL){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  d <- d/param
  return((1+sqrt(3)*d)*exp(-sqrt(3)*d))
}

matern5_hsic <- function(x,param,d=NULL){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  d <- d/param
  return((1+sqrt(5)*d+5/3*d^2)*exp(-sqrt(5)*d))
}

# Bernoulli polynomials
B1 <- function(t){
  return(t-.5)
}
B2 <- function(t){
  return(t^2-t+1/6)
}
B4 <- function(t){
  return(t^4-2*t^3+t^2-1/30)
}

ssanova1_hsic <- function(x,d=NULL,...){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  return(1+B1(x)%*%t(B1(x))+0.5*B2(d))
}

ssanova2_hsic <- function(x,d=NULL,...){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  return(1+B1(x)%*%t(B1(x))+B2(x)%*%t(B2(x))/4-B4(d)/24)
}

# Old version for computing the estimators - 2014
# HSIC <- function(X, Y, kernelX, paramX, kernelY, paramY, estimator.type = "V-stat"){
#   
#   nobs <- nrow(X)
#   
#   Centering matrix
#   H <- diag(nobs) - 1/nobs*matrix(1,nobs,nobs)
# 
#   KX <- do.call(get(paste(kernelX,"_hsic",sep="")), list(x=X,param=paramX))
#   KXH <- H%*%KX%*%H
# 
#   KY <- do.call(get(paste(kernelY,"_hsic",sep="")), list(x=Y,param=paramY))
#   KYH <- H%*%KY%*%H
# 
#   return(list(estimate=mean(KXH*KYH),meanX=mean(KXH*KXH),meanY=mean(KYH*KYH)))
# }

# ---------------------------------------------------------------------------------------------- 
# New HSIC function - 2019
HSIC <- function(X, Y, kernelX, paramX, kernelY, paramY, estimator.type = "V-stat"){
  # Two types of estimators available : 
  # if estimator.type == "V-stat" (default value): biased (but asymptotically unbiased) estimator, more practical for numerical implementation
  # if estimator.type == "U-stat": unbiased estimator
  # The variance is of order o(1/n) for both estimators
  # For more details, see Meynaoui et al. (2019)
  
  nobs <- nrow(X)  #sample size
  d <- ncol(X)     #number of scalar inputs
  HSICXY <- rep(0,d)
  HSICXX <- rep(0,d)
  
  KY <- do.call(get(paste(kernelY,"_hsic",sep="")), list(x=Y,param=paramY))
  
  if(estimator.type == "V-stat"){
    
    UYY <- sum(KY*KY)
    VYY <- sum(colSums(KY)%*%KY)
    WYY <- sum(KY)*sum(KY)
    HSICYY <- UYY/nobs^2 - 2*VYY/nobs^3 + WYY/nobs^4
    
    for (r in 1:d) {
      KX <- do.call(get(paste(kernelX,"_hsic",sep="")), list(x=X[,r],param=paramX)) 
      UXY <- sum(KX*KY)
      UXX <- sum(KX*KX)
      
      VXY <- sum(colSums(KX)%*%KY)
      VXX <- sum(colSums(KX)%*%KX)
      
      WXY <- sum(KX)*sum(KY)
      WXX <- sum(KX)*sum(KX)
      
      HSICXY[r] <- UXY/nobs^2 - 2*VXY/nobs^3 + WXY/nobs^4
      HSICXX[r] <- UXX/nobs^2 - 2*VXX/nobs^3 + WXX/nobs^4
    }
  } 
  
  if(estimator.type == "U-stat"){
    
    diag(KY) <- 0 #Put all the Y-diag elements to zero
    
    UYY <- sum(KY*KY)
    VYY <- sum(colSums(KY)%*%KY) - UYY
    WYY <- sum(KY)*sum(KY) - 2*UYY - 4*VYY
  
    HSICYY <- UYY/(nobs*(nobs-1)) - 2*VYY/(nobs*(nobs-1)*(nobs-2)) + WYY/(nobs*(nobs-1)*(nobs-2)*(nobs-3))
    
    for (r in 1:d) {
      KX <- do.call(get(paste(kernelX,"_hsic",sep="")), list(x=X[,r],param=paramX)) 
      diag(KX) <- 0 #Put all the X-diag elements to zero
      
      UXY <- sum(KX*KY)
      UXX <- sum(KX*KX)
      
      VXY <- sum(colSums(KX)%*%KY) - UXY
      VXX <- sum(colSums(KX)%*%KX) - UXX
      
      WXY <- sum(KX)*sum(KY) - 2*UXY - 4*VXY
      WXX <- sum(KX)*sum(KX) - 2*UXX - 4*VXX
      
      HSICXY[r] <- UXY/(nobs*(nobs-1)) - 2*VXY/(nobs*(nobs-1)*(nobs-2)) + WXY/(nobs*(nobs-1)*(nobs-2)*(nobs-3))
      HSICXX[r] <- UXX/(nobs*(nobs-1)) - 2*VXX/(nobs*(nobs-1)*(nobs-2)) + WXX/(nobs*(nobs-1)*(nobs-2)*(nobs-3))
    } 
  }
  return(list(estimate=HSICXY,meanX=HSICXX,meanY=HSICYY))
}

# ---------------------------------------------------------------------------------------------- 
sensiHSIC <- function(model = NULL, X, kernelX = "rbf", paramX = NA, 
                      kernelY = "rbf", paramY = NA, nboot = 0, conf = 0.95,
                      estimator.type = "V-stat", test.method = "Asymptotic", B = 1000, ...) {
  
  
  if (is.data.frame(X)){
    X <- as.matrix(unname(X))
  }
  else if(!is.matrix(X)){
    stop("The sample X must be a matrix or a data frame")
  }
  
  p <- ncol(X)
  nkx <- length(kernelX)
  if (!(nkx == 1 | nkx == p)){
    stop("KernelX must be of length 1 or p (number of input variables)")
  }
  if (!missing(paramX)){
    npx <- length(paramX)
    if (!(npx == 1 | npx == p)){
      stop("paramX must be of length 1 or p (number of input variables)")
    }
  }
  
  if(test.method == "Asymptotic"){
    x <- list(model = model, X = X, kernelX = kernelX, paramX = paramX, 
              kernelY = kernelY, paramY = paramY, nboot = nboot,
              conf = conf, estimator.type = estimator.type, test.method = test.method , call = match.call()) 
  }else{
    x <- list(model = model, X = X, kernelX = kernelX, paramX = paramX, 
              kernelY = kernelY, paramY = paramY, nboot = nboot,
              conf = conf, estimator.type = estimator.type, test.method = estimator.type, B = B, call = match.call()) 
  }
  class(x) <- "sensiHSIC"
  
  #calcul of the response for explicit model
  if (! is.null(x$model)) {
    response(x, ...)
    x=tell(x, ...)
  }  
  return(x)
}


estim.sensiHSIC <- function(data, i=1:nrow(data), kernelX, paramX, 
                            kernelY, paramY, estimator.type) {
  
  ptot <- ncol(data)
  p <- ptot - 1
  X <- data[i,1:p]
  Y <- data[i,ptot]
  S = matrix(0,nrow=p,ncol=1)
  HSICXY = matrix(0,nrow=p,ncol=1)
  
  # HSIC indices
  for (i in 1:p){
    Xtemp <- as.matrix(X[,i])
    Ytemp <- as.matrix(Y)
    res <- HSIC(Xtemp, Ytemp, kernelX[i], paramX[i], kernelY, paramY, estimator.type = estimator.type)
    
    HSICXY[i] <- res$estimate   # HSIC
    S[i] <- res$estimate/sqrt(res$meanX*res$meanY)  # R2 HSIC
    
  }
  return(list(S=S, HSICXY = HSICXY))
}

# ----------------------------------------------------------------------------------------------  
# Function to compute a p-value of independence test based on HSIC measure as statistic
estim.sensiHSIC.pvalue <- function(data, kernelX, paramX, 
                                   kernelY, paramY, estimator.type,test.method, B) {
  
  ptot <- ncol(data)
  p <- ptot - 1
  X <- data[,1:p]
  Y <- data[,ptot]
  P = matrix(0,nrow=p,ncol=1)
  
  # HSIC indices
  for (i in 1:p){
    Xtemp <- as.matrix(X[,i])
    Ytemp <- as.matrix(Y)
    
    if(test.method == "Permutation"){#P-value computed by permutation-based method with B bootstrap samples
      P[i] <- perm_test_HSIC(Xtemp, Ytemp, kernelX[i], paramX[i], kernelY, paramY,estimator.type = estimator.type, B = B)
    }
    else{# P-value computed by asymtotic approximation (default method)
      P[i] <- asymp_test_HSIC(Xtemp, Ytemp, kernelX[i], paramX[i], kernelY, paramY,estimator.type =estimator.type)
    }
  }
  return(P)
}
# ---------------------------------------------------------------------------------------------- 
# Two elementary functions to access to each output of estim.sensiHSIC function (required to compute the confidence interval with bootstrap)
estim.sensiHSIC.S <- function(data,i=1:nrow(data), kernelX, paramX, 
                              kernelY, paramY, estimator.type) {
  res = estim.sensiHSIC(data=data, i= i, kernelX, paramX,kernelY, paramY,estimator.type = estimator.type)$S
  return(res)
}

estim.sensiHSIC.HSICXY <- function(data, i=1:nrow(data), kernelX, paramX, 
                                   kernelY, paramY, estimator.type) {
  res = estim.sensiHSIC(data=data, i= i, kernelX, paramX,kernelY, paramY,estimator.type = estimator.type)$HSICXY
  return(res)
}
# ---------------------------------------------------------------------------------------------- 

tell.sensiHSIC <- function(x, y = NULL, ...) {
  
  id <- deparse(substitute(x))
  if (! is.null(y)) {
    x$y <- y
  } 
  else if (is.null(x$y)) {
    stop("y not found")
  }
  
  n <- nrow(x$X)
  p <- ncol(x$X)
  
  nkx <- length(x$kernelX)
  if (nkx==1) x$kernelX <- rep(x$kernelX,p)
  
  if (is.na(x$paramX[1])){
    x$paramX <- matrix(nrow=1,ncol=p)
    for (i in 1:p){
      if (x$kernelX[i]=="dcov"){
        x$paramX[i] = 1
      }
      else{
        x$paramX[i] = sd(x$X[,i])
      }
    }
  }else{
    if (length(x$paramX)==1) x$paramX <- rep(x$paramX,p)
  }
  
  
  if (is.na(x$paramY)){
    if (x$kernelY=="dcov"){
      x$paramY = 1
    }else{
      x$paramY <- sd(x$y)
    }
  }
  
  data <-cbind(x$X,x$y)
  
  if (x$nboot == 0) {
    res <- estim.sensiHSIC(data,1:n, kernelX = x$kernelX, paramX = x$paramX,  
                           kernelY = x$kernelY, paramY =x$paramY, estimator.type = x$estimator.type)
    x$S <- data.frame(res$S)
    colnames(x$S) <- "original"
    x$HSICXY <- data.frame(res$HSICXY)
    colnames(x$HSICXY) <- "original"
    
  } 
  else {
    S.boot <- boot(data, estim.sensiHSIC.S, 
                   kernelX = x$kernelX, paramX = x$paramX, 
                   kernelY = x$kernelY, paramY = x$paramY, 
                   estimator.type = x$estimator.type, R = x$nboot)
    
    x$S <- bootstats(S.boot, x$conf, "basic")
    
    HSICXY.boot <- boot(data, estim.sensiHSIC.HSICXY, 
                        kernelX = x$kernelX, paramX = x$paramX, 
                        kernelY = x$kernelY, paramY = x$paramY, estimator.type = x$estimator.type, R = x$nboot)
    x$HSICXY <- bootstats(HSICXY.boot, x$conf, "basic")
  }
  
  # Computation of pvalue from independence tests with H0 : independence hypothesis
  res <- estim.sensiHSIC.pvalue(data, x$kernelX, x$paramX, x$kernelY, x$paramY, estimator.type = x$estimator.type,test.method = x$test.method, B = x$B)
  x$Pvalue <- data.frame(res)
  colnames(x$Pvalue) <- "original"
  
  rownames <- paste("X",1:p,sep="")
  rownames(x$S) <- rownames(x$HSICXY) <- rownames(x$Pvalue) <- rownames 
  
  assign(id, x, parent.frame())
  return(x)
}


print.sensiHSIC<- function(x, ...) {
  
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (! is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    if (! is.null(x$S)) {
      cat("\n\n\nSensitivity HSIC indices (denoted S or R2HSIC) \n")
      print(x$S)
      if (! is.null(x$HSICXY)) { cat("\n\n\nHSIC indices\n")
        print(x$HSICXY)
      }
      if (! is.null(x$Pvalue)) { cat("\n\n\nP-value\n")
        print(x$Pvalue)
      }
    }
    else{
      cat("(empty)\n")
    }
  }
}


plot.sensiHSIC <- function(x, ylim = c(0, 1), ...) {
  
  if (! is.null(x$y)) {
    nodeplot(x$S, ylim = ylim)
    legend(x = "topright", legend = "HSIC Sensitivity Indices")
  }
}

ggplot.sensiHSIC <- function(x, ylim = c(0, 1), ...) {
  
  if (! is.null(x$y)) {
    nodeggplot(list(x$S), xname = "HSIC Sensitivity indices", ylim = ylim)
  }
}

# --------------------------------------------------------------------------
# Functions to compute the pvalue of independence test, with HISC measure as statistic test.
# H0: X and Y are independent (<=> HSIC (X,Y) = 0 with universal kernels) against H1: X and Y are dependent 
# The P-value (under H0) can be computed:
#        - for asymptotic framework: with asymptotic approximation (Gamma approximation) => asymp_test_HSIC function
#        - for non-asymptotic framewok: with permutation approach, based on B boostrap samples => perm_test_HSIC function
# For more details: Meynaoui et al. (2019) - New statistical methodology for second-level global sensitivity analysis of numerical simulators.
#
#    Anouar Meynaoui & Amandine Marrel 2019
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# Asymptotic test of independence with HISC measure as statistic test.
# X: Input matrix or data.frame, number of columns = number of scalar inputs,  number of rows = sample size 
# Y: Output vector or matrix or data frame, number of columns = number of scalar outputs,  number of rows = sample size
# kernelX: the chosen kernel for the inputs
# paramX: the parameters of kernelX
# kernelY: the chosen kernel for the outputs
# paramY: the parameters of kernelY
# estimator.type : the type of estimation, "V-stat" for V-statistics or "U-stat" for U-statistics, see HSIC function and Meynaoui et al. (2019) for details

asymp_test_HSIC <- function(X, Y, kernelX = "rbf", paramX = sd(X), kernelY = "rbf", paramY = sd(Y), estimator.type  = "V-stat"){
  
  if (is.data.frame(X) | is.vector(X)) {
    X <- as.matrix(unname(X))
  }
  else if (!is.matrix(X)) {
    stop("The sample X must be a vector or a matrix or a data frame")
  }
  
  d <- ncol(X) #number of scalar inputs
  n <- nrow(X) #sample size
  
  if (is.data.frame(Y) | is.vector(Y)) {
    Y <- as.matrix(unname(Y))
  }
  else if (!is.matrix(Y)) {
    stop("The sample Y must be a vector or matrix or a data frame")
  }
  
  if (nrow(Y) != n) {
    stop("The sample size of the inputs and the outputs is not the same")
  }
  
  pval <- rep(0,d) #The vector of pvalues to be calculated by this function
  
  H <- diag(n)-matrix(1/n,n,n) #Useful matrix from [Gretton2008] to estimate the shape and scale of the Gamma distribution 
  KY <- do.call(get(paste(kernelY,"_hsic",sep="")), list(x=Y,param=paramY)) #Kernel matrix of Y
  
  muy <- 1/(n*(n-1))*sum(KY-diag(KY)*diag(ncol(KY)))
  BY <- H%*%KY%*%H
  
  for(s in 1:d){
    KX <- do.call(get(paste(kernelX,"_hsic",sep="")), list(x=X[,s],param=paramX)) #kernel matrix of the input X_s
    
    mux <- 1/(n*(n-1))*sum(KX-diag(KX)*diag(ncol(KX)))
    BX <- H%*%KX%*%H
    
    B <- (BX*BY)^2 
    
    HSMean <- 1/n*(1 + mux*muy - mux - muy)
    HSVariance <- (2*(n-4)*(n-5)/(n*(n-1)*(n-2)*(n-3)))*sum(B-diag(B)*diag(ncol(B)))/n/(n-1)
    
    alpha <- (HSMean)^2/HSVariance
    beta  <- n*HSVariance/(HSMean)
    
    if(estimator.type  == "U-stat"){# U-stat estimator (unbiased)
      nHSIC <- n*HSIC(as.matrix(X[,s]),cbind(Y),kernelX, paramX, kernelY, paramY, estimator.type  = "U-stat")$estimate
      pval[s] <- 1-pgamma(q = nHSIC + n*HSMean , shape = alpha, rate = 1/beta)
    }else{# V-stat estimator by default (biased bu asymptotically unbiased and more practical for implementation)
      nHSIC <- n*HSIC(as.matrix(X[,s]),cbind(Y),kernelX, paramX, kernelY, paramY, estimator.type  = "V-stat")$estimate
      pval[s] <- 1-pgamma(q = nHSIC, shape = alpha, rate = 1/beta)}
  }
  return(pval)
}

# --------------------------------------------------------------------------
# Permutation-based test of independence with HISC measure as statistic test.
# X: Input matrix or data.frame, number of columns = number of scalar inputs,  number of rows = sample size 
# Y: Output vector or matrix or data frame, number of columns = number of scalar outputs,  number of rows = sample size
# kernelX: the chosen kernel for the inputs
# paramX: the parameters of kernelX
# kernelY: the chosen kernel for the outputs
# paramY: the parameters of kernelY
# estimator.type : the type of estimation, "V-stat" for V-statistics or "U-stat" for U-statistics, see HSIC function and Meynaoui et al. (2019) for details
# B: number of permutations by default B = 1000

perm_test_HSIC <- function(X, Y, kernelX = "rbf", paramX = sd(X), kernelY = "rbf", paramY = sd(Y), estimator.type  = "V-stat",B = 1000){
  
  if (is.data.frame(X) | is.vector(X)) {
    X <- as.matrix(unname(X))
  }
  else if (!is.matrix(X)) {
    stop("The sample X must be a vector or a matrix or a data frame")
  }
  
  d <- ncol(X) #number of scalar inputs
  n <- nrow(X) #sample size
  
  if (is.data.frame(Y) | is.vector(Y)) {
    Y <- as.matrix(unname(Y))
  }
  else if (!is.matrix(Y)) {
    stop("The sample Y must be a vector or matrix or a data frame")
  }
  
  if (nrow(Y) != n) {
    stop("The sample size of the inputs and the outputs is not the same")
  }
  
  pval <- rep(0,d) #The vector of pvalues to be calculated by this function
  
  KY <- do.call(get(paste(kernelY,"_hsic",sep="")), list(x=Y,param=paramY)) #Kernel matrix of Y
  
  if(estimator.type  == "U-stat"){  
    
    Estimated.value <- HSIC(X, Y, kernelX, paramX, kernelY, paramY, estimator.type  = "U-stat")$estimate
    
    diag(KY) <- 0 #Put all the diagonal elements to zero
    
    for(s in 1:d){
      KX <- do.call(get(paste(kernelX,"_hsic",sep="")), list(x=X[,s],param=paramX)) #kernel matrix of the input X_s
      diag(KX) <- 0 #Put all the diagonal elements to zero
      
      
      hsic.perm <- sapply(1:B,function(al){ #permutations under the null hypothesis of HSIC estimator
        shuf <- sample(n) 
        KYY <- KY[shuf,shuf]
        U <- sum(KYY*KX) #We use this element twice
        V <- sum(colSums(KX)%*%KYY)-U #We use this element twice
        W <- sum(KYY)*sum(KX) - 4*V -2*U
        U/(n*(n-1))+ W/(n*(n-1)*(n-2)*(n-3))-2*V/(n*(n-1)*(n-2))})
      
      hsic.perm <- c(hsic.perm,Estimated.value[s]) #The new estimation method from [Meynaoui et al., 2019]  
      
      pval[s] <- length(hsic.perm[hsic.perm > Estimated.value[s]])/length(hsic.perm)} 
  }
  else{  # V-stat estimator by default
    Estimated.value <- HSIC(X, Y, kernelX, paramX, kernelY, paramY, estimator.type  = "V-stat")$estimate
    
    for(s in 1:d){
      KX <- do.call(get(paste(kernelX,"_hsic",sep="")), list(x=X[,s],param=paramX)) #kernel matrix of the input X_s
      
      hsic.perm <- sapply(1:B,function(al){ #permutations under the null hypothesis of HSIC estimator
        shuf <- sample(n) 
        KYY <- KY[shuf,shuf]
        U <- sum(KYY*KX) 
        V <- sum(colSums(KX)%*%KYY) 
        W <- sum(KYY)*sum(KX)
        U/n^2 + W/n^4 - 2*V/n^3})
      
      hsic.perm <- c(hsic.perm,Estimated.value[s]) #The new estimation method from [Meynaoui et al., 2019]  
      
      pval[s] <- length(hsic.perm[hsic.perm > Estimated.value[s]])/length(hsic.perm)} 
  }
  
  return(pval)
}


