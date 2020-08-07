##################################################################
# Sensitivity indices based on kernel embeddings of distributions
#
# Sebastien Da Veiga 2014, Anouar Meynaoui & Amandine Marrel 2019
##################################################################

my_fun <- function(a,b) {
  if (!requireNamespace("pkg", quietly = TRUE)) {
    stop("Package \"pkg\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
}

# Kernel functions
rbf_hsic <- function(x,param,d=NULL){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  d <- d/param
  return(exp(-0.5*d^2))
}

categ_hsic <- function(x,param,d=NULL){
  if (is.null(d)){
    if (param==0){
      d <- 1*sapply(x,function(z){z==x})
    }else{
      n <- nrow(x)
      d <- matrix(0,n,n)
      val <- unique(x)
      n.val <- length(val)
      for(i in 1:n.val){
        id <- which(x==val[i])
        d[id,id] <- 1/sum(x==val[i])
      }
    }
  }
  return(d)
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

  q <- ncol(as.matrix(Y))     #number of scalar outputs
  KY <- 1
  for (i in 1:q){
    KY <- KY * do.call(get(paste(kernelY,"_hsic",sep="")), list(x=Y[,i,drop=FALSE],param=paramY[i]))
  }

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
# New C-HSIC function - 2020
CHSIC <- function(X, Y, cond, kernelX, paramX, kernelY, paramY, estimator.type = "V-stat"){
  # Available estimator:
  # if estimator.type == "V-stat" (default value): biased (but asymptotically unbiased) estimator, more practical for numerical implementation
  # For more details, see Marrel (2020)

  trace <- function(data) sum(diag(data))
  nobs <- nrow(X)  #sample size
  d <- ncol(X)     #number of scalar inputs
  HSICXY <- rep(0,d)
  HSICXX <- rep(0,d)
  if(!is.null(cond)){ #weight calculation
    weight <- weightTSA(Y, cond$c, cond$upper, cond$type, cond$param)
    weights <- diag(weight/mean(weight))
  }
  KY <- do.call(get(paste(kernelY,"_hsic",sep="")), list(x=Y,param=paramY))

  if(estimator.type == "V-stat"){

    #H1 <- diag(nobs) - 1/nobs*(matrix(1,nobs,nobs)%*%weights)
    #H2 <- diag(nobs) - 1/nobs*(weights%*%matrix(1,nobs,nobs))
    #HSICYY <- trace(weights%*%KY%*%weights%*%H1%*%KY%*%H2)/(nobs^2)
    w1 <- t(apply(weights,1,function(x) rep(sum(x),length(x))))
    w2 <- (apply(weights,1,function(x) rep(sum(x),length(x))))
    W <- w1*w2
    UYY <- sum(KY*KY*W)
    VYY <- sum((colSums(KY*w1)%*%(KY*W)))
    WYY <- sum(KY*W)*sum(KY*W)
    HSICYY <- UYY/nobs^2 - 2*VYY/nobs^3 + WYY/nobs^4

    for (r in 1:d) {
      KX <- do.call(get(paste(kernelX,"_hsic",sep="")), list(x=X[,r],param=paramX))

      #HSICXY[r] <- trace(weights%*%KX%*%weights%*%H1%*%KY%*%H2)/(nobs^2)
      UXY <- sum(KX*KY*W)
      VXY <- sum((colSums(KX*w1)%*%(KY*W)))
      WXY <- sum(KX*W)*sum(KY*W)
      HSICXY[r] <- UXY/nobs^2 - 2*VXY/nobs^3 + WXY/nobs^4

      #HSICXX[r] <- trace(weights%*%KX%*%weights%*%H1%*%KX%*%H2)/(nobs^2)
      UXX <- sum(KX*KX*W)
      VXX <- sum((colSums(KX*w1)%*%(KX*W)))
      WXX <- sum(KX*W)*sum(KX*W)
      HSICXX[r] <- UXX/nobs^2 - 2*VXX/nobs^3 + WXX/nobs^4
    }
  }
  else{
    # w1 <- t(apply(weights,1,function(x) rep(sum(x),length(x))))
    # w2 <- (apply(weights,1,function(x) rep(sum(x),length(x))))
    # W <- w1*w2
    # diag(KY) <- 0 #Put all the Y-diag elements to zero
    # UYY <- sum(KY*KY*W)
    # VYY <- sum(colSums(KY*w1)%*%(KY*W)) - UYY
    # WYY <- sum(KY*W)*sum(KY*W) - 2*UYY - 4*VYY
    # HSICYY <- UYY/(nobs*(nobs-1)) + WYY/(nobs*(nobs-1)*(nobs-2)*(nobs-3)) - 2*VYY/(nobs*(nobs-1)*(nobs-2))

    nobsQ <- length(which(weight!=0))
    ID <- which(weight!=0)
    KY_Q <- do.call(get(paste(kernelY = "rbf","_hsic",sep="")), list(x=Y[ID],param=paramY))
    diag(KY_Q) <- 0 #Put all the Y-diag elements to zero
    UYY <- sum(KY_Q*KY_Q)
    VYY <- sum(colSums(KY_Q)%*%(KY_Q)) - UYY
    WYY <- sum(KY_Q)*sum(KY_Q) - 2*UYY - 4*VYY
    HSICYY <- UYY/(nobsQ*(nobsQ-1)) - 2*VYY/(nobsQ*(nobsQ-1)*(nobsQ-2)) +
      WYY/(nobsQ*(nobsQ-1)*(nobsQ-2)*(nobsQ-3))

    for (r in 1:d) {
      # KX <- do.call(get(paste(kernelX,"_hsic",sep="")), list(x=X[,r],param=paramX))
      # diag(KX) <- 0 #Put all the X-diag elements to zero
      # UXX <- sum(KX*KX*W)
      # VXX <- sum(colSums(KX*w1)%*%(KX*W)) - UXX
      # WXX <- sum(KX*W)*sum(KX*W) - 2*UXX - 4*VXX
      # HSICXX[r] <- UXX/(nobs*(nobs-1)) + WXX/(nobs*(nobs-1)*(nobs-2)*(nobs-3)) - 2*VXX/(nobs*(nobs-1)*(nobs-2))
      # UXY <- sum(KX*KY*W)
      # VXY <- sum(colSums(KX*w1)%*%(KY*W)) - UXY
      # WXY <- sum(KX*W)*sum(KY*W) - 2*UXY - 4*VXY
      # HSICXY[r] <- UXY/(nobs*(nobs-1)) + WXY/(nobs*(nobs-1)*(nobs-2)*(nobs-3)) - 2*VXY/(nobs*(nobs-1)*(nobs-2))

      KX_Q <- do.call(get(paste(kernelX = "rbf","_hsic",sep="")), list(x=X[ID,r],param=paramX))
      diag(KX_Q) <- 0 #Put all the X-diag elements to zero
      UXY <- sum(KX_Q*KY_Q) #sum(KX*KY*W)
      UXX <- sum(KX_Q*KX_Q) #sum(KX*KX*W)
      VXY <- sum(colSums(KX_Q)%*%(KY_Q)) - UXY #sum(colSums(KX*w1)%*%(KY*W)) - UXY
      VXX <- sum(colSums(KX_Q)%*%(KX_Q)) - UXX #sum(colSums(KX*w1)%*%(KX*W)) - UXX
      WXY <- sum(KX_Q)*sum(KY_Q) - 2*UXY - 4*VXY #sum(KX*W)*sum(KY*W) - 2*UXY - 4*VXY
      WXX <- sum(KX_Q)*sum(KX_Q) - 2*UXX - 4*VXX #sum(KX*W)*sum(KX*W) - 2*UXX - 4*VXX
      HSICXY[r] <- UXY/(nobsQ*(nobsQ-1)) - 2*VXY/(nobsQ*(nobsQ-1)*(nobsQ-2)) +
        WXY/(nobsQ*(nobsQ-1)*(nobsQ-2)*(nobsQ-3))
      HSICXX[r] <- UXX/(nobsQ*(nobsQ-1)) - 2*VXX/(nobsQ*(nobsQ-1)*(nobsQ-2)) +
        WXX/(nobsQ*(nobsQ-1)*(nobsQ-2)*(nobsQ-3))
    }
  }
  return(list(estimate=HSICXY,meanX=HSICXX,meanY=HSICYY))
}

# ----------------------------------------------------------------------------------------------
sensiHSIC <- function(model = NULL, X, target = NULL, cond = NULL,
                      kernelX = "rbf", paramX = NA,
                      kernelY = "rbf", paramY = NA, nboot = 0, conf = 0.95,
                      estimator.type = "V-stat", test.method = "Asymptotic", B = 5000,
                      crit.option = list(stop.criterion = "screening", alpha = 0.05, Bstart = 100,
                                         Bfinal = 5000, Bbatch = 100, lo = 200, graph = TRUE), expl.var.PCA = NULL, ...) {
  if(!(estimator.type == "V-stat" | estimator.type == "U-stat")){
    estimator.type = "V-stat"
    warning("estimator.type must be V-stat or U-stat. By default V-stat has been performed")
  }
  if(!(test.method == "Asymptotic" | test.method == "Permutation" | test.method == "Seq_Permutation" | test.method == "No")){
    test.method = "Permutation"
    warning("test.method must be Asymptotic, Permutation, Seq_Permutation or No. By default Permutation has been performed")
  }
  if(!is.null(cond) && (test.method == "Asymptotic")){
    warning("Asymptotic test is not yet implemented. By default Permutation has been performed.")
    test.method = "Permutation"
  }
  if(!is.null(target) && target$type == "indicTh" && kernelY != "categ"){
    warning("A more appropriate kernel should be used for binary output. By default, categorical kernel has been applied.")
    kernelY = "categ"
  }
  if(is.null(crit.option$stop.criterion)) crit.option$stop.criterion <- "screening"
  if(is.null(crit.option$alpha)) crit.option$alpha <- 0.05
  if(is.null(crit.option$Bstart)) crit.option$Bstart <- 100
  if(is.null(crit.option$Bfinal)) crit.option$Bfinal <- 5000
  if(is.null(crit.option$Bbatch)) crit.option$Bbatch <- 100
  if(is.null(crit.option$lo)) crit.option$lo <- 200
  if(is.null(crit.option$graph)) crit.option$graph <- TRUE

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
    x <- list(model = model, X = X, target = target, cond = cond, kernelX = kernelX, paramX = paramX,
              kernelY = kernelY, paramY = paramY, nboot = nboot,
              conf = conf, estimator.type = estimator.type, test.method = test.method ,
              expl.var.PCA = expl.var.PCA, call = match.call())
  }
  else if(test.method == "No"){
    x <- list(model = model, X = X, target = target, cond = cond, kernelX = kernelX, paramX = paramX,
              kernelY = kernelY, paramY = paramY, nboot = nboot,
              conf = conf, estimator.type = estimator.type, test.method = test.method ,
              expl.var.PCA = expl.var.PCA, call = match.call())
  }
  else if(test.method == "Seq_Permutation"){
    x <- list(model = model, X = X, target = target, cond = cond, kernelX = kernelX, paramX = paramX,
              kernelY = kernelY, paramY = paramY, nboot = nboot,
              conf = conf, estimator.type = estimator.type, test.method = test.method,
              crit.option = crit.option, expl.var.PCA = expl.var.PCA, call = match.call())
  }
  else{
    x <- list(model = model, X = X, target = target, cond = cond, kernelX = kernelX, paramX = paramX,
              kernelY = kernelY, paramY = paramY, nboot = nboot,
              conf = conf, estimator.type = estimator.type, test.method = test.method,
              B = B, expl.var.PCA = expl.var.PCA, call = match.call())
  }
  class(x) <- "sensiHSIC"

  #calcul of the response for explicit model
  if (!is.null(x$model)) {
    response(x, ...)
    if(!is.null(target)){
      if (is.null(target$c)) stop("threshold not found")
      if(is.null(target$type)) target$type <- "indicTh"
      if(is.null(target$upper)) target$upper <- TRUE
      if(target$type == "zeroTh"){
        target$type = "indicTh"
        warning("Weight function must be indicTh, logistic or exp1side. By default indicTh has been performed")
      }
      if(is.null(target$param)) target$param <- 1
      x$y <- weightTSA(x$y, target$c, target$upper, target$type, target$param)
    }
    if(!is.null(cond)){
      if(is.null(cond$c)) stop("threshold not found")
      if(is.null(cond$type)) cond$type <- "exp1side"
      if(is.null(cond$upper)) cond$upper <- TRUE
      if(is.null(cond$param)) cond$param <- 1
    }
    x=tell(x, ...)
  }
  return(x)
}

estim.sensiHSIC <- function(data, i=1:nrow(data), cond, kernelX, paramX,
                            kernelY, paramY, estimator.type, q) {
  ptot <- ncol(data)
  p <- ptot - q
  X <- as.matrix(data[i,1:p],ncol=p)
  Y <- data[i,(p+1):ptot,drop=FALSE]
  S = matrix(0,nrow=p,ncol=1)
  HSICXY = matrix(0,nrow=p,ncol=1)
  # HSIC indices
  for (i in 1:p){
    Xtemp <- as.matrix(X[,i])
    Ytemp <- as.matrix(Y)
    if(is.null(cond))
      res <- HSIC(Xtemp, Ytemp, kernelX[i], paramX[i], kernelY, paramY, estimator.type = estimator.type)
    else
      res <- CHSIC(Xtemp, Ytemp, cond, kernelX[i], paramX[i], kernelY, paramY, estimator.type = estimator.type)

    HSICXY[i] <- res$estimate   # HSIC
    S[i] <- res$estimate/sqrt(res$meanX*res$meanY)  # R2 HSIC
  }
  return(list(S=S, HSICXY = HSICXY))
}

# ----------------------------------------------------------------------------------------------
# Function to compute a p-value of independence test based on HSIC measure as statistic
estim.sensiHSIC.pvalue <- function(data, kernelX, paramX,
                                   kernelY, paramY, estimator.type,test.method, B, crit.option, q) {
  ptot <- ncol(data)
  p <- ptot - q
  X <- as.matrix(data[,1:p],ncol=p)
  Y <- data[,(p+1):ptot,drop=FALSE]
  P = matrix(0,nrow=p,ncol=1)
  PP <- NULL
  if(p==1 && test.method == "Seq_Permutation"){
    if(crit.option$stop.criterion == "ranking" | crit.option$stop.criterion == "both"){
      crit.option$stop.criterion = "screening"
      warning("Only screening can be used in the case of univariate input variable. By default screening has been performed")
    }
  }

  if(test.method == "Permutation"){
    for (i in 1:p){
      Xtemp <- as.matrix(X[,i])
      Ytemp <- as.matrix(Y)
      P[i] <- perm_test_HSIC(Xtemp, Ytemp, kernelX[i], paramX[i], kernelY, paramY,
                             estimator.type = estimator.type, B = B)$pval
    }
  }

  else if(test.method == "Seq_Permutation"){
    Estimated.value <- NULL
    for(i in 1:p)
      Estimated.value[i] <- HSIC(as.matrix(X[,i]), as.matrix(Y), kernelX[i], paramX[i], kernelY, paramY,
                                 estimator.type  = estimator.type)$estimate
    if(!(crit.option$stop.criterion == "screening" | crit.option$stop.criterion == "ranking" | crit.option$stop.criterion == "both")){
      crit.option$stop.criterion = "screening"
      warning("Type of stop.criterion must be screening, ranking or both. By default screening has been performed")
    }
    if(crit.option$Bfinal < crit.option$lo) stop("lo must be lower than Bfinal")
    if(crit.option$alpha < 0 || crit.option$alpha > 1) stop("alpha should be in the interval [0,1]")
    if(crit.option$Bstart > crit.option$Bfinal) stop("Bstart must lower than Bfinal")
    if(crit.option$Bstart + crit.option$Bbatch > crit.option$Bfinal){
      for (i in 1:p){
        Xtemp <- as.matrix(X[,i])
        Ytemp <- as.matrix(Y)
        P[i] <- perm_test_HSIC(Xtemp, Ytemp, kernelX[i], paramX[i], kernelY, paramY,
                               estimator.type = estimator.type, B = crit.option$Bfinal)$pval
      }
      #P <- formatC(P, format = "e", digits = 5)
    }
    if(crit.option$Bstart + crit.option$Bbatch <= crit.option$Bfinal  && crit.option$stop.criterion != "screening")
    {
      B_seq <- seq(crit.option$Bstart+crit.option$Bbatch,crit.option$Bfinal,by=crit.option$Bbatch)
      seq_HSIC <- matrix(NA,nrow=p,ncol=crit.option$Bfinal)

      for(i in 1:p){
        Xtemp <- as.matrix(X[,i])
        Ytemp <- as.matrix(Y)
        p_eq <- perm_test_HSIC(Xtemp, Ytemp, kernelX[i], paramX[i], kernelY, paramY,estimator.type = estimator.type, B = crit.option$Bstart)$hsic.perm[-(crit.option$Bstart+1)]
        seq_HSIC[i,1:length(p_eq)] <- p_eq
      }

      for(k in 1:length(B_seq))
      {
        for(i in 1:p)
        {
          Xtemp <- as.matrix(X[,i])
          Ytemp <- as.matrix(Y)
          p_eq <- perm_test_HSIC(Xtemp, Ytemp, kernelX[i], paramX[i], kernelY, paramY,estimator.type = estimator.type, B = crit.option$Bbatch)$hsic.perm[-(crit.option$Bbatch+1)]
          muet <- seq_HSIC[i,]
          seq_HSIC[i,1:length(c(muet[!is.na(muet)],p_eq))] <- c(muet[!is.na(muet)],p_eq)
        }
        pvalue_sequence <- matrix(NA,nrow=p,ncol=B_seq[k]-1)
        for(m in 1:(B_seq[k]-1)) for(i in 1:p) {val <- c(seq_HSIC[i,1:(m+1)],Estimated.value[i]) ;  pvalue_sequence[i,m] <- length(val[val > Estimated.value[i]])/length(val) }

        pvalue_screening <- pvalue_sequence
        pvalue_screening[pvalue_screening <= crit.option$alpha] <- 0 ; pvalue_screening[pvalue_screening > crit.option$alpha] <- 1

        if(crit.option$stop.criterion == "ranking")
          if(dim(pvalue_sequence)[2] >= crit.option$lo &&
             all(apply(apply(pvalue_sequence[,(dim(pvalue_sequence)[2]- crit.option$lo):dim(pvalue_sequence)[2]],2,rank,ties.method ="first"),1,diff) == 0) == TRUE) break

        if(crit.option$stop.criterion == "both")
          if(dim(pvalue_sequence)[2] >= crit.option$lo &&
             all(apply(pvalue_screening[,(dim(pvalue_screening)[2]- crit.option$lo):dim(pvalue_screening)[2]],1,diff) == 0) == TRUE &&
             all(apply(apply(pvalue_sequence[,(dim(pvalue_sequence)[2]- crit.option$lo):dim(pvalue_sequence)[2]],2,rank,ties.method ="first"),1,diff) == 0) == TRUE) break
      }
      PP <- t(apply(pvalue_sequence,1,as.numeric)) #pvalue_sequence
      P <- pvalue_sequence[1:p,dim(pvalue_sequence)[2]]
    }
    if(crit.option$Bstart + crit.option$Bbatch <= crit.option$Bfinal && crit.option$stop.criterion == "screening")
    {
      B_seq <- seq(crit.option$Bstart+crit.option$Bbatch,crit.option$Bfinal,by=crit.option$Bbatch)
      seq_HSIC <- matrix(NA,nrow=p,ncol=crit.option$Bfinal)
      PP <- matrix(NA,nrow=p,ncol=crit.option$Bfinal)

      for(i in 1:p){
        Xtemp <- as.matrix(X[,i])
        Ytemp <- as.matrix(Y)
        p_eq <- perm_test_HSIC(Xtemp, Ytemp, kernelX[i], paramX[i], kernelY, paramY,estimator.type = estimator.type, B = crit.option$Bstart)$hsic.perm[-(crit.option$Bstart+1)]
        seq_HSIC[i,1:length(p_eq)] <- p_eq
      }

      for(i in 1:p)
      {
        seq_HSIC_1D <- matrix(seq_HSIC[i,],nrow=1)
        for(k in 1:length(B_seq))
        {
          Xtemp <- as.matrix(X[,i])
          Ytemp <- as.matrix(Y)
          p_eq <- perm_test_HSIC(Xtemp, Ytemp, kernelX[i], paramX[i], kernelY, paramY,
                                 estimator.type = estimator.type, B = crit.option$Bbatch)$hsic.perm[-(crit.option$Bbatch+1)]
          muet <- seq_HSIC_1D[1,]
          seq_HSIC_1D[1,1:length(c(muet[!is.na(muet)],p_eq))] <- c(muet[!is.na(muet)],p_eq)

          pvalue_sequence <- matrix(NA,nrow=1,ncol=B_seq[k]-1)
          dim(pvalue_sequence)
          for(m in 1:(B_seq[k]-1)){
            val <- c(seq_HSIC_1D[1,1:(m+1)],Estimated.value[i])
            pvalue_sequence[1,m] <- length(val[val > Estimated.value[i]])/length(val)
          }

          pvalue_screening <- matrix(pvalue_sequence,nrow=1)
          pvalue_screening[pvalue_screening <= 0.05] <- 0
          pvalue_screening[pvalue_screening > 0.05] <- 1
          if( dim(pvalue_screening)[2] >= crit.option$lo ){
            if(all(as.numeric(pvalue_screening[1,(dim(pvalue_screening)[2]- crit.option$lo):dim(pvalue_screening)[2]]) == rep(1,crit.option$lo+1)) |
               all(as.numeric(pvalue_screening[1,(dim(pvalue_screening)[2]- crit.option$lo):dim(pvalue_screening)[2]]) == rep(0,crit.option$lo+1))) break
          }
        }
        P[i] <- pvalue_sequence[length(pvalue_sequence)]
        PP[i,1:length(pvalue_sequence)] <- pvalue_sequence
      }
    }
  }

  else{
    for (i in 1:p){
      Xtemp <- as.matrix(X[,i])
      Ytemp <- as.matrix(Y)
      P[i] <- asymp_test_HSIC(Xtemp, Ytemp, kernelX[i], paramX[i], kernelY, paramY,estimator.type = estimator.type)
    }
  }
  return(list(P=P,Pseq=PP))
}

# ----------------------------------------------------------------------------------------------
# Function to compute a p-value of independence test based on C-HSIC measure as statistic
estim.sensiCHSIC.pvalue <- function(data, cond, kernelX, paramX,
                                    kernelY, paramY, estimator.type,test.method, B, crit.option) {
  ptot <- ncol(data)
  p <- ptot - 1
  X <- as.matrix(data[,1:p],ncol=p)
  Y <- data[,ptot]
  P = matrix(0,nrow=p,ncol=1)
  PP <- NULL

  if(p==1 && test.method == "Seq_Permutation"){
    if(crit.option$stop.criterion == "ranking" | crit.option$stop.criterion == "both"){
      crit.option$stop.criterion = "screening"
      warning("Only screening can be used in the case of univariate input variable. By default screening has been performed")
    }
  }

  if(test.method == "Permutation"){
    for (i in 1:p){
      Xtemp <- as.matrix(X[,i])
      Ytemp <- as.matrix(Y)
      P[i] <- perm_test_CHSIC(Xtemp, Ytemp, cond, kernelX[i], paramX[i], kernelY, paramY,
                              estimator.type = estimator.type, B = B)$pval
    }
  }

  else if(test.method == "Seq_Permutation"){
    Estimated.value <- NULL
    for(i in 1:p)
      Estimated.value[i] <- CHSIC(as.matrix(X[,i]), as.matrix(Y), cond, kernelX[i], paramX[i], kernelY, paramY,
                                  estimator.type  = estimator.type)$estimate

    if(!(crit.option$stop.criterion == "screening" | crit.option$stop.criterion == "ranking" | crit.option$stop.criterion == "both")){
      crit.option$stop.criterion = "screening"
      warning("Type of stop.criterion must be screening or ranking. By default screening has been performed")
    }
    if(crit.option$Bfinal < crit.option$lo) stop("lo must be lower than Bfinal")
    if(crit.option$alpha < 0 || crit.option$alpha > 1) stop("alpha must be lower than 1")
    if(crit.option$Bstart > crit.option$Bfinal) stop("Bstart must lower than Bfinal")
    if(crit.option$Bstart + crit.option$Bbatch > crit.option$Bfinal){
      for (i in 1:p){
        Xtemp <- as.matrix(X[,i])
        Ytemp <- as.matrix(Y)
        P[i] <- perm_test_CHSIC(Xtemp, Ytemp, cond, kernelX[i], paramX[i], kernelY, paramY,
                                estimator.type = estimator.type, B = crit.option$Bfinal)$pval
      }
      #P <- formatC(P, format = "e", digits = 5)
    }
    if(crit.option$Bstart + crit.option$Bbatch <= crit.option$Bfinal  && crit.option$stop.criterion != "screening")
    {
      B_seq <- seq(crit.option$Bstart+crit.option$Bbatch,crit.option$Bfinal,by=crit.option$Bbatch)
      seq_HSIC <- matrix(NA,nrow=p,ncol=crit.option$Bfinal)

      for(i in 1:p){
        Xtemp <- as.matrix(X[,i])
        Ytemp <- as.matrix(Y)
        p_eq <- perm_test_CHSIC(Xtemp, Ytemp, cond, kernelX[i], paramX[i], kernelY, paramY,estimator.type = estimator.type, B = crit.option$Bstart)$hsic.perm[-(crit.option$Bstart+1)]
        seq_HSIC[i,1:length(p_eq)] <- p_eq
      }

      for(k in 1:length(B_seq))
      {
        for(i in 1:p)
        {
          Xtemp <- as.matrix(X[,i])
          Ytemp <- as.matrix(Y)
          p_eq <- perm_test_CHSIC(Xtemp, Ytemp, cond, kernelX[i], paramX[i], kernelY, paramY,estimator.type = estimator.type, B = crit.option$Bbatch)$hsic.perm[-(crit.option$Bbatch+1)]
          muet <- seq_HSIC[i,]
          seq_HSIC[i,1:length(c(muet[!is.na(muet)],p_eq))] <- c(muet[!is.na(muet)],p_eq)
        }
        pvalue_sequence <- matrix(NA,nrow=p,ncol=B_seq[k]-1)
        for(m in 1:(B_seq[k]-1)) for(i in 1:p) {val <- c(seq_HSIC[i,1:(m+1)],Estimated.value[i]) ;  pvalue_sequence[i,m] <- length(val[val > Estimated.value[i]])/length(val) }

        pvalue_screening <- pvalue_sequence
        pvalue_screening[pvalue_screening <= crit.option$alpha] <- 0 ; pvalue_screening[pvalue_screening > crit.option$alpha] <- 1

        if(crit.option$stop.criterion == "ranking")
          if(dim(pvalue_sequence)[2] >= crit.option$lo &&
             all(apply(apply(pvalue_sequence[,(dim(pvalue_sequence)[2]- crit.option$lo):dim(pvalue_sequence)[2]],2,rank,ties.method ="first"),1,diff) == 0) == TRUE) break

        if(crit.option$stop.criterion == "both")
          if(dim(pvalue_sequence)[2] >= crit.option$lo &&
             all(apply(pvalue_screening[,(dim(pvalue_screening)[2]- crit.option$lo):dim(pvalue_screening)[2]],1,diff) == 0) == TRUE &&
             all(apply(apply(pvalue_sequence[,(dim(pvalue_sequence)[2]- crit.option$lo):dim(pvalue_sequence)[2]],2,rank,ties.method ="first"),1,diff) == 0) == TRUE) break
      }
      PP <- t(apply(pvalue_sequence,1,as.numeric)) #pvalue_sequence
      P <- pvalue_sequence[1:p,dim(pvalue_sequence)[2]]
    }
    if(crit.option$Bstart + crit.option$Bbatch <= crit.option$Bfinal && crit.option$stop.criterion == "screening")
    {
      B_seq <- seq(crit.option$Bstart+crit.option$Bbatch,crit.option$Bfinal,by=crit.option$Bbatch)
      seq_HSIC <- matrix(NA,nrow=p,ncol=crit.option$Bfinal)
      PP <- matrix(NA,nrow=p,ncol=crit.option$Bfinal)

      for(i in 1:p){
        Xtemp <- as.matrix(X[,i])
        Ytemp <- as.matrix(Y)
        p_eq <- perm_test_CHSIC(Xtemp, Ytemp, cond, kernelX[i], paramX[i], kernelY, paramY,estimator.type = estimator.type, B = crit.option$Bstart)$hsic.perm[-(crit.option$Bstart+1)]
        seq_HSIC[i,1:length(p_eq)] <- p_eq
      }

      for(i in 1:p)
      {
        seq_HSIC_1D <- matrix(seq_HSIC[i,],nrow=1)
        for(k in 1:length(B_seq))
        {
          Xtemp <- as.matrix(X[,i])
          Ytemp <- as.matrix(Y)
          p_eq <- perm_test_CHSIC(Xtemp, Ytemp, cond, kernelX[i], paramX[i], kernelY, paramY,
                                  estimator.type = estimator.type, B = crit.option$Bbatch)$hsic.perm[-(crit.option$Bbatch+1)]
          muet <- seq_HSIC_1D[1,]
          seq_HSIC_1D[1,1:length(c(muet[!is.na(muet)],p_eq))] <- c(muet[!is.na(muet)],p_eq)

          pvalue_sequence <- matrix(NA,nrow=1,ncol=B_seq[k]-1)
          dim(pvalue_sequence)
          for(m in 1:(B_seq[k]-1)){
            val <- c(seq_HSIC_1D[1,1:(m+1)],Estimated.value[i])
            pvalue_sequence[1,m] <- length(val[val > Estimated.value[i]])/length(val)
          }

          pvalue_screening <- matrix(pvalue_sequence,nrow=1)
          pvalue_screening[pvalue_screening <= 0.05] <- 0
          pvalue_screening[pvalue_screening > 0.05] <- 1

          if( dim(pvalue_screening)[2] >= crit.option$lo ){
            print(all(as.numeric(pvalue_screening[1,(dim(pvalue_screening)[2]- crit.option$lo):dim(pvalue_screening)[2]]) == rep(1,crit.option$lo+1)))
            print(all(as.numeric(pvalue_screening[1,(dim(pvalue_screening)[2]- crit.option$lo):dim(pvalue_screening)[2]]) == rep(0,crit.option$lo+1)))
            if(all(as.numeric(pvalue_screening[1,(dim(pvalue_screening)[2]- crit.option$lo):dim(pvalue_screening)[2]]) == rep(1,crit.option$lo+1)) |
               all(as.numeric(pvalue_screening[1,(dim(pvalue_screening)[2]- crit.option$lo):dim(pvalue_screening)[2]]) == rep(0,crit.option$lo+1))) break
          }
        }
        P[i] <- pvalue_sequence[length(pvalue_sequence)]
        PP[i,1:length(pvalue_sequence)] <- pvalue_sequence
      }
    }
  }

  else{
    for (i in 1:p){
      Xtemp <- as.matrix(X[,i])
      Ytemp <- as.matrix(Y)
      P[i] <- asymp_test_CHSIC(Xtemp, Ytemp, cond, kernelX[i], paramX[i], kernelY, paramY,estimator.type = estimator.type)
    }
  }
  return(list(P=P,Pseq=PP))
}

# ----------------------------------------------------------------------------------------------
# Two elementary functions to access to each output of estim.sensiHSIC function (required to compute the confidence interval with bootstrap)
estim.sensiHSIC.S <- function(data,i=1:nrow(data), cond, kernelX, paramX,
                              kernelY, paramY, estimator.type, q) {
  res = estim.sensiHSIC(data=data, i= i, cond, kernelX, paramX,kernelY, paramY,estimator.type = estimator.type, q = q)$S
  return(res)
}

estim.sensiHSIC.HSICXY <- function(data, i=1:nrow(data), cond, kernelX, paramX,
                                   kernelY, paramY, estimator.type, q) {
  res = estim.sensiHSIC(data=data, i= i, cond, kernelX, paramX,kernelY, paramY,estimator.type = estimator.type, q = q)$HSICXY
  return(res)
}

# ----------------------------------------------------------------------------------------------

tell.sensiHSIC <- function(x, y = NULL, ...) {
  id <- deparse(substitute(x))
  if (!is.null(y) && is.null(x$target)) {
    x$y <- y
  }
  if (!is.null(y) && !is.null(x$target)){
    if (is.null(x$target$c)) stop("threshold not found")
    if(is.null(x$target$type)) x$target$type <- "indicTh"
    if(is.null(x$target$upper)) x$target$upper <- TRUE
    if(x$target$type == "zeroTh"){
      x$target$type = "indicTh"
      warning("Weight function must be indicTh, logistic or exp1side. By default indicTh has been performed")
    }
    if(is.null(x$target$param)) x$target$param <- 1
    x$y <- weightTSA(y, x$target$c, x$target$upper, x$target$type, x$target$param)
  }
  else if (is.null(x$y)) {
    stop("y not found")
  }
  n <- nrow(x$X)
  p <- ncol(x$X)

  if (is.null(dim(x$y))){
    x$y <- matrix(x$y,ncol=1)
  }else{
    x$y <- as.matrix(unname(x$y))
  }
  q <- ncol(x$y)

  if (q>1 & any(x$kernelY=="categ")){
    stop("Cannot use a categorical kernel with multiple outputs yet")
  }

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

  if (q>1){ # Multiple outputs
    if (!is.null(x$expl.var.PCA)){ # We use a PCA transform on the outputs
      # First remove constant outputs, if any
      id.cst <- apply(x$y,2,sd,na.rm=TRUE)==0
      # Perform PCA
      pca.res <- prcomp(x$y[,!id.cst], retx = TRUE, center = TRUE, scale. = TRUE)
      # Keep components until explained variance is greater than expl.var.PCA
      expl.var <- cumsum(pca.res$sdev^2)/sum(pca.res$sdev^2)
      k.pca <- min(which(expl.var>x$expl.var.PCA))
      x$y <- unname(pca.res$x[,1:k.pca])
      q <- k.pca
    }
    nky <- length(x$kernelY)
    if (nky==1) x$kernelY <- rep(x$kernelY,q)
    if (is.na(x$paramY[1])){
      x$paramY <- matrix(nrow=1,ncol=q)
      for (i in 1:q){
        if (x$kernelY[i]=="dcov"){
          x$paramY[i] = 1
        }
        else{
          x$paramY[i] = sd(x$y[,i])
        }
      }
    }else{
      if (length(x$paramY)==1) x$paramY <- rep(x$paramY,q)
    }
  }else{
    # Categorical kernel is treated differently
    if (x$kernelY=="categ"){
      y.val <- unique(x$y)
      ny.val <- length(y.val)
      if (is.na(x$paramY[1]) | x$paramY[1]=="normal"){
        x$paramY <- 0
      }else if (x$paramY[1]=="weighted"){
        x$paramY <- 1
      }else{
        stop("paramY for categorical kernel is not equal to normal or weighted")
      }
    }else if (is.na(x$paramY[1])){
      if (x$kernelY=="dcov"){
        x$paramY = 1
      }else{
        x$paramY = sd(x$y)
      }
    }
  }

  data <- cbind(x$X,x$y)

  if (x$nboot == 0) {
    res <- estim.sensiHSIC(data, 1:n, cond = x$cond, kernelX = x$kernelX, paramX = x$paramX,
                           kernelY = x$kernelY, paramY =x$paramY, estimator.type = x$estimator.type, q = q)
    x$S <- data.frame(res$S)
    colnames(x$S) <- "original"
    x$HSICXY <- data.frame(res$HSICXY)
    colnames(x$HSICXY) <- "original"
  }
  else {
    S.boot <- boot(data, estim.sensiHSIC.S,
                   cond = x$cond, kernelX = x$kernelX, paramX = x$paramX,
                   kernelY = x$kernelY, paramY = x$paramY,
                   estimator.type = x$estimator.type, q = q, R = x$nboot)
    x$S <- bootstats(S.boot, x$conf, "basic")
    HSICXY.boot <- boot(data, estim.sensiHSIC.HSICXY,
                        cond = x$cond, kernelX = x$kernelX, paramX = x$paramX,
                        kernelY = x$kernelY, paramY = x$paramY, estimator.type = x$estimator.type, q = q, R = x$nboot)
    x$HSICXY <- bootstats(HSICXY.boot, x$conf, "basic")
  }

  # Computation of pvalue from independence tests with H0 : independence hypothesis

  if(x$test.method != "No"){
    if(is.null(x$cond)){
      res <- estim.sensiHSIC.pvalue(data, x$kernelX, x$paramX, x$kernelY, x$paramY,
                                    estimator.type = x$estimator.type,test.method = x$test.method,
                                    B = x$B, crit.option = x$crit.option, q = q)
    }
    else{
      res <- estim.sensiCHSIC.pvalue(data, x$cond, x$kernelX, x$paramX, x$kernelY, x$paramY,
                                     estimator.type = x$estimator.type,test.method = x$test.method,
                                     B = x$B, crit.option = x$crit.option)
    }
  }


  if(x$test.method == "Seq_Permutation"){
    x$Pvalue <- data.frame((res$P))
    len <- apply(res$Pseq,1,function(x) length(na.omit(x)))
    x$SeqPvalue <- data.frame((res$Pseq[,(1:max(len))]))
    colnames(x$Pvalue) <- "original"
    rownames <- paste("X",1:p,sep="")
    rownames(x$S) <- rownames(x$HSICXY) <- rownames(x$Pvalue) <- rownames
    colnames(x$SeqPvalue) <- NULL
    assign(id, x, parent.frame())
  }
  if(x$test.method == "Seq_Permutation" && x$crit.option$graph == TRUE)
    matplot(t(x$SeqPvalue),type='l',xlab="Number of permutations",
            ylab="The Pvalues",main="Sequential Estimation of the Pvalues",ylim=c(0,1))
  if(x$test.method == "No"){
    rownames <- paste("X",1:p,sep="")
    rownames(x$S) <- rownames(x$HSICXY) <- rownames
    assign(id, x, parent.frame())
  }
  else{
    x$Pvalue <- data.frame(res$P)
    colnames(x$Pvalue) <- "original"
    rownames <- paste("X",1:p,sep="")
    rownames(x$S) <- rownames(x$HSICXY) <- rownames(x$Pvalue) <- rownames
    assign(id, x, parent.frame())
  }
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

plot.sensiHSIC <- function(x, ylim = c(0,1), ...) {

  if (! is.null(x$y)) {
    nodeplot(x$S, ylim = ylim)
    legend(x = "topright", legend = "HSIC Sensitivity Indices")
  }
}

ggplot.sensiHSIC <- function(x, ylim = c(0,1), ...) {

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

  q <- ncol(as.matrix(Y))     #number of scalar outputs
  KY <- 1
  for (i in 1:q){
    KY <- KY * do.call(get(paste(kernelY,"_hsic",sep="")), list(x=Y[,i,drop=FALSE],param=paramY[i]))
  }

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

perm_test_HSIC <- function(X, Y, kernelX = "rbf", paramX = sd(X), kernelY = "rbf", paramY = sd(Y), estimator.type  = "V-stat",B = B){

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

  q <- ncol(as.matrix(Y))     #number of scalar outputs
  KY <- 1
  for (i in 1:q){
    KY <- KY * do.call(get(paste(kernelY,"_hsic",sep="")), list(x=Y[,i,drop=FALSE],param=paramY[i]))
  }

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
        shuf <- sample(seq(1,n),n) #sample(n)
        KYY <- KY[shuf,shuf]
        U <- sum(KYY*KX)
        V <- sum(colSums(KX)%*%KYY)
        W <- sum(KYY)*sum(KX)
        U/n^2 + W/n^4 - 2*V/n^3
      }
      )
      hsic.perm <- c(hsic.perm,Estimated.value[s]) #The new estimation method from [Meynaoui et al., 2019]
      pval[s] <- length(hsic.perm[hsic.perm > Estimated.value[s]])/length(hsic.perm)}
  }

  return(list(pval=pval,hsic.perm=hsic.perm))
}

# --------------------------------------------------------------------------
# Asymptotic test of independence with Cond-HISC measure as statistic test.
# X: Input matrix or data.frame, number of columns = number of scalar inputs,  number of rows = sample size
# Y: Output vector or matrix or data frame, number of columns = number of scalar outputs,  number of rows = sample size
# kernelX: the chosen kernel for the inputs
# paramX: the parameters of kernelX
# kernelY: the chosen kernel for the outputs
# paramY: the parameters of kernelY
# estimator.type : the type of estimation, "V-stat" for V-statistics or "U-stat" for U-statistics, see HSIC function and Meynaoui et al. (2019) for details

asymp_test_CHSIC <- function(X, Y, cond, kernelX = "rbf", paramX = sd(X), kernelY = "rbf", paramY = sd(Y), estimator.type  = "V-stat"){

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

  weight <- weightTSA(Y, cond$c, cond$upper, cond$type, cond$param)
  weights <- diag(weight/mean(weight))

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
# Permutation-based test of independence with Cond-HISC measure as statistic test.
# X: Input matrix or data.frame, number of columns = number of scalar inputs,  number of rows = sample size
# Y: Output vector or matrix or data frame, number of columns = number of scalar outputs,  number of rows = sample size
# cond: list of parameters to perform Conditional Sensitivity Analysis (CSA)
# kernelX: the chosen kernel for the inputs
# paramX: the parameters of kernelX
# kernelY: the chosen kernel for the outputs
# paramY: the parameters of kernelY
# estimator.type : the type of estimation, "V-stat" for V-statistics or "U-stat" for U-statistics, see HSIC function and Meynaoui et al. (2019) for details
# B: number of permutations by default B = 1000

perm_test_CHSIC <- function(X, Y, cond, kernelX = "rbf", paramX = sd(X), kernelY = "rbf", paramY = sd(Y),
                            estimator.type  = "V-stat",B = B){

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

  weight <- weightTSA(Y, cond$c, cond$upper, cond$type, cond$param)
  weights <- diag(weight/mean(weight))
  w1 <- t(apply(weights,1,function(x) rep(sum(x),length(x))))
  w2 <- (apply(weights,1,function(x) rep(sum(x),length(x))))
  W <- w1*w2
  nobsQ <- (length(which(weight!=0)))

  if(estimator.type  == "U-stat"){

    Estimated.value <- CHSIC(X, Y, cond, kernelX, paramX, kernelY, paramY, estimator.type  = "U-stat")$estimate

    for(s in 1:d){

      hsic.perm <- sapply(1:B,function(al){ #permutations under the null hypothesis of HSIC estimator

        shuf <- sample(seq(1,n),n) #sample(n)
        #CHSIC(X, Y[shuf], cond, kernelX, paramX, kernelY, paramY, estimator.type  = "U-stat")$estimate
        Y_shuf <- Y[shuf]
        weight_shuf <- weightTSA(Y_shuf, cond$c, cond$upper, cond$type, cond$param)
        ID_shuf <- which(weight_shuf!=0)
        nobsQ <- (length(which(weight_shuf!=0)))
        KY_Q <- do.call(get(paste(kernelY = "rbf","_hsic",sep="")), list(x=Y_shuf[ID_shuf],param=paramY))
        diag(KY_Q) <- 0 #Put all the Y-diag elements to zero
        KX_Q <- do.call(get(paste(kernelX = "rbf","_hsic",sep="")), list(x=X[ID_shuf,s],param=paramX))
        diag(KX_Q) <- 0 #Put all the X-diag elements to zero
        UXY <- sum(KX_Q*KY_Q)
        VXY <- sum(colSums(KX_Q)%*%(KY_Q)) - UXY
        WXY <- sum(KX_Q)*sum(KY_Q) - 2*UXY - 4*VXY
        UXY/(nobsQ*(nobsQ-1)) - 2*VXY/(nobsQ*(nobsQ-1)*(nobsQ-2)) + WXY/(nobsQ*(nobsQ-1)*(nobsQ-2)*(nobsQ-3))
      })
      hsic.perm <- c(hsic.perm,Estimated.value[s])
      pval[s] <- length(hsic.perm[hsic.perm > Estimated.value[s]])/length(hsic.perm)}
  }
  else{

    KY <- do.call(get(paste(kernelY,"_hsic",sep="")), list(x=Y,param=paramY)) #Kernel matrix of Y
    Estimated.value <- CHSIC(X,Y, cond, kernelX, paramX, kernelY, paramY, estimator.type  = "V-stat")$estimate

    for(s in 1:d){
      KX <- do.call(get(paste(kernelX,"_hsic",sep="")), list(x=X[,s],param=paramX)) #kernel matrix of the input X_s
      hsic.perm <- sapply(1:B,function(al){ #permutations under the null hypothesis of HSIC estimator
        shuf <- sample(seq(1,n),n) #sample(n)
        WW <- W[shuf,shuf]
        ww1 <- w1[shuf,shuf]
        KYY <- KY[shuf,shuf]
        UXY <- sum(KX*KYY*WW)
        VXY <- sum((colSums(KX*ww1)%*%(KYY*WW)))
        WXY <- sum(KX*WW)*sum(KYY*WW)
        UXY/n^2 + WXY/n^4 - 2*VXY/n^3})

      hsic.perm <- c(hsic.perm,Estimated.value[s]) #The new estimation method from [Meynaoui et al., 2019]
      pval[s] <- length(hsic.perm[hsic.perm > Estimated.value[s]])/length(hsic.perm)}
  }
  return(list(pval=pval,hsic.perm=hsic.perm))
}

