### testHSIC ###################################################################

testHSIC <- function(sensi, 
                     test.method="Asymptotic",
                     B=3000, seq.options=list(criterion="screening", alpha=0.05, 
                                              Bstart=200, Bfinal=5000, Bbatch=100, Bconv=200, 
                                              graph=TRUE) 
                     ){
  
  ### inputs ###

  # sensi:            1 x ...     list of class sensiHSIC (with all data)
  # test.method:      1 x 1       character (method used to test independence)
  # B:                1 x 1       numeric   (number of permutations)
  # seq.options:      1 x ...     list      (options for the sequential test)

  ### output ###

  # x:                1 x ...     list of class testHSIC (with all results)

  ### useful lists #############################################################
  
  L <- available.material()

  ### For sensi ################################################################

  if(missing(sensi)){
    stop("You must specificy what sensi is.", call.=FALSE)
  }

  ### For test.method ##########################################################

  # test.method must belong to L$tests

  if(!is.character(test.method) | length(test.method)>1){
    stop("test.method must be a string.", call.=FALSE)
  }else{
    if(!(test.method %in% L$tests)){
      stop("Invalid test method. \n",
           "Available test methods include: ",
           paste0(L$tests, collapse=', '), ".", call.=FALSE)
    }
  }
  
  # two statistics can be used to test independence:
  # --> the numerator of the first-order HSIC-ANOVA index
  # --> the numerator of the total-order HSIC-ANOVA index
  
  if(stringr::str_detect(test.method, "^Tot_")){
    measure.def <- "TO.num" 
  }else{
    measure.def <- "HSICXY"
  }
  
  ### For sensi and test.method ################################################

  check.sensiHSIC(sensi, test.method)

  ### For B ####################################################################

  if(!is.numeric(B) | length(B)>1){
    stop("B must be a positive integer.", call.=FALSE)
  }else{
    if(B<=0){
      stop("B must be a positive integer.", call.=FALSE)
    }else{
      if(!(B%%1==0)){
        B <- ceiling(B)
        warning("B was not an integer. \n",
                "ceiling(B) has been taken instead.", call.=FALSE)
      }
    }
  }

  ### For seq.options ##########################################################

  ### minimum requirement ###

  # 1) seq.options must be a list of options
  # 2) seq.options must contain at least one element named "criterion"

  if(!is.list(seq.options)){
    stop("seq.options must be a list of options.", call.=FALSE)
  }else{
    if(is.null(seq.options$criterion)){
      stop("Criterion not defined. You must specify seq.options$criterion.", call.=FALSE)
    }
  }

  ### default values for other options ###

  if(is.null(seq.options$alpha)) seq.options$alpha <- 0.05
  if(is.null(seq.options$Bstart)) seq.options$Bstart <- 200
  if(is.null(seq.options$Bfinal)) seq.options$Bfinal <- 5000
  if(is.null(seq.options$Bbatch)) seq.options$Bbatch <- 100
  if(is.null(seq.options$Bconv)) seq.options$Bconv <- 200
  if(is.null(seq.options$graph)) seq.options$graph <- TRUE

  ### removal of undesired options in seq.options ###

  expected.seq <- c("criterion", "alpha",
                    "Bstart", "Bfinal", "Bbatch", "Bconv", "graph")

  found.seq <- names(seq.options)
  if(!(all(found.seq %in% expected.seq))){
    undesired.seq <- setdiff(found.seq, expected.seq)
    seq.options[undesired.seq] <- NULL
    warning("The list seq.options can only contain the following options: \n",
            paste0(expected.seq, collapse=', '), ". \n",
            "Other options have been removed.", call.=FALSE)
  }

  ### check of remaining options in seq.options ###

  # seq.options$criterion must be a string

  if(!is.character(seq.options$criterion) | length(seq.options$criterion)>1){
    stop("seq.options$criterion must be a string.", call.=FALSE)
  }else{
    if(!(seq.options$criterion %in% c("screening", "ranking", "both"))){
      seq.options$criterion <- "screening"
      warning("seq.options$criterion must be: screening, ranking or both. \n",
              "By default, seq.options$criterion=screening has been selected.", call.=FALSE)
    }
  }

  # seq.options$alpha must be a real number in [0,1]

  if(!is.numeric(seq.options$alpha) | length(seq.options$alpha)>1){
    stop("seq.options$alpha must be a real number between 0 and 1.", call.=FALSE)
  }else{
    if(seq.options$alpha<0 | seq.options$alpha>1){
      seq.options$alpha <- 0.05
      warning("seq.options$alpha was not between 0 and 1. \n",
              "By default, seq.options$alpha=0.05 has been selected.", call.=FALSE)
    }
  }

  # seq.options$graph must be boolean

  if(length(seq.options$graph)>1){
    stop("seq.options$graph must be a boolean.", call.=FALSE)
  }else{
    if(!(isTRUE(seq.options$graph) | isFALSE(seq.options$graph))){
      stop("seq.options$graph must be a boolean.", call.=FALSE)
    }
  }

  ### special focus on Bstart, Bfinal, Bbatch and Bconv ###

  # the hyperparameters are unpacked
  
  Bstart <- seq.options$Bstart
  Bfinal <- seq.options$Bfinal
  Bbatch <- seq.options$Bbatch
  Bconv <- seq.options$Bconv

  # Bstart must be a positive integer

  if(!is.numeric(Bstart) | length(Bstart)>1){
    stop("seq.options$Bstart must be a positive integer.", call.=FALSE)
  }else{
    if(Bstart<=0){
      stop("seq.options$Bstart must be a positive integer.", call.=FALSE)
    }else{
      if(!(Bstart%%1==0)){
        Bstart <- ceiling(Bstart)
        warning("seq.options$Bstart was not an integer. \n",
                "ceiling(seq.options$Bstart) has been taken instead.", call.=FALSE)
      }
    }
  }

  # Bfinal must be a positive integer

  if(!is.numeric(Bfinal) | length(Bfinal)>1){
    stop("seq.options$Bfinal must be a positive integer.", call.=FALSE)
  }else{
    if(Bfinal<=0){
      stop("seq.options$Bfinal must be a positive integer.", call.=FALSE)
    }else{
      if(!(Bfinal%%1==0)){
        Bfinal <- ceiling(Bfinal)
        warning("seq.options$Bfinal was not an integer. \n",
                "ceiling(seq.options$Bfinal) has been taken instead.", call.=FALSE)
      }
    }
  }

  # Bbatch must be a positive integer

  if(!is.numeric(Bbatch) | length(Bbatch)>1){
    stop("seq.options$Bbatch must be a positive integer.", call.=FALSE)
  }else{
    if(Bbatch<=0){
      stop("seq.options$Bbatch must be a positive integer.", call.=FALSE)
    }else{
      if(!(Bbatch%%1==0)){
        Bbatch <- ceiling(Bbatch)
        warning("seq.options$Bbatch was not an integer. \n",
                "ceiling(seq.options$Bbatch) has been taken instead.", call.=FALSE)
      }
    }
  }

  # Bconv must be a positive integer

  if(!is.numeric(Bconv) | length(Bconv)>1){
    stop("seq.options$Bconv must be a positive integer.", call.=FALSE)
  }else{
    if(Bconv<=0){
      stop("seq.options$Bconv must be a positive integer.", call.=FALSE)
    }else{
      if(!(Bconv%%1==0)){
        Bconv <- ceiling(Bconv)
        warning("seq.options$Bconv was not an integer. \n",
                "ceiling(seq.options$Bconv) has been taken instead.", call.=FALSE)
      }
    }
  }

  ### additional precautions ###
  
  need.replacement <- FALSE
  
  # precaution 1: Bstart must be greater than Bfinal
  
  if(Bstart>=Bfinal){
    stop("Bfinal must be larger than Bstart. \n",
         "Here, you chose Bstart=", Bstart, 
         " and Binal=", Bfinal, ".", call.=FALSE)
  }

  # precaution 2: Bconv must be greater than or equal to Bfinal
  
  if(Bconv>Bfinal){
    stop("Bfinal must be larger than or equal to Bconv. \n",
         "Here, you chose Bconv=", Bconv, 
         " and Binal=", Bfinal, ".", call.=FALSE)
  }

  # precaution 3: at least one batch may be added
  
  if(Bbatch>Bfinal-Bstart){
    Bbatch <- Bfinal-Bstart
    need.replacement <- TRUE
  }

  # precaution 4: Bfinal must be equal to Bstart+k*Bbatch

  if(!((Bfinal-Bstart)%%Bbatch==0)){
    cat("boolean=",!((Bfinal-Bstart)%%Bbatch==0))
    Bfinal <- Bstart+ceiling((Bfinal-Bstart)/Bbatch)*Bbatch
    need.replacement <- TRUE
  }

  # to send a warning if some values have been corrected
  
  if(need.replacement){
    
    warning("Some hyperparameters in seq.options are not consistent. \n",
            "The following values have been selected: \n",
            "Bstart=", Bstart, ", ", "Bfinal=", Bfinal, ", ",
            "Bbatch=", Bbatch, ", ", "Bconv=", Bconv, ".", call.=FALSE)
    
    # the hyperparameters are repacked
    
    seq.options$Bstart <- Bstart
    seq.options$Bfinal <- Bfinal
    seq.options$Bbatch <- Bbatch
    seq.options$Bconv <- Bconv

  }
  
  #########################
  ### Main instructions ###
  #########################
  
  ### extraction from sensi ###
  
  HSIC.obs <- switch(measure.def, 
                     "HSICXY"={ unname(sensi$HSICXY[,"original"]) },
                     "TO.num"={ unname(sensi$TO.num[,"original"]) },
                     NA)
  
  ### specific case of C-HSIC indices ###
  
  if(!is.null(sensi$cond)){
    mask.cond <- sensi$weights>0
    KX <- sensi$KX[mask.cond, mask.cond, ]
    KY <- sensi$KY[mask.cond, mask.cond]
  }else{
    KX <- sensi$KX
    KY <- sensi$KY
  }
  
  ### independence test ###

  if(test.method=="Asymptotic"){
    
    res <- asymptotic.test(HSIC.obs, KX, KY, sensi$estimator.type)
    
  }else if(test.method %in% c("Gamma", "Tot_Gamma")){
    
    res <- gamma.test(HSIC.obs, KX, KY, measure.def)
    
  }else if(test.method %in% c("Permutation", "Tot_Permutation")){
    
    res <- permutation.test(HSIC.obs, KX, KY, 
                            measure.def, sensi$estimator.type, B) 
    
  }else if(test.method %in% c("Seq_Permutation", "Tot_Seq_Permutation")){
    
    res <- seq.permutation.test(HSIC.obs, KX, KY,
                                measure.def, sensi$estimator.type, seq.options)
    
  }
  
  ### object of class testHSIC ###
  
  x <- list()
  
  p <- dim(KX)[3]
  
  # p-values
  x$pval <- data.frame(res$pval)
  rownames(x$pval) <- paste("X", 1:p, sep="")
  colnames(x$pval) <- "original"
  
  # properties of the chosen test
  
  x$prop <- switch(test.method,
                   "Asymptotic"={ c("asymptotic", "parametric") }, 
                   "Permutation"={ c("non-asymptotic", "non-parametric") },
                   "Seq_Permutation"={ c("non-asymptotic", "non-parametric") },
                   "Gamma"={ c("non-asymptotic", "parametric") },
                   "Tot_Permutation"={ c("non-asymptotic", "non-parametric") },
                   "Tot_Seq_Permutation"={ c("non-asymptotic", "non-parametric") },
                   "Tot_Gamma"={ c("non-asymptotic", "parametric") },
                   NA)
  
  # parametric family
  
  if(x$prop[2]=="parametric"){
    
    if(test.method=="Asymptotic"){
      x$family <- switch(sensi$estimator.type, 
                         "U-stat"={ "Pearson3" }, 
                         "V-stat"={ "Gamma" }, NA)
    }else{
      x$family <- "Gamma"
    }
    
  }
  
  # about parameters (for parametric methods)
  
  if(x$prop[2]=="parametric"){
  
    x$param <- data.frame(t(res$param))
    rownames(x$param) <- paste("X", 1:p, sep="")
    
    if(test.method=="Asymptotic"){
      colnames(x$param) <- switch(sensi$estimator.type,
                                  "U-stat"={ c("shape", "scale", "position") }, 
                                  "V-stat"={ c("shape", "scale") }, NA)
    }else{
      colnames(x$param) <- c("shape", "scale")
    }
    
  }

  # about permutations (for non-parametric methods)
  
  if(test.method %in% c("Permutation", "Tot_Permutation")){
    
    x$B <- B
    
    x$Hperm <- data.frame(t(res$Hperm))
    rownames(x$Hperm) <- paste("X", 1:p, sep="")
    colnames(x$Hperm) <- NULL
    
  }
  
  if(test.method %in% c("Seq_Permutation", "Tot_Seq_Permutation")){
    
    x$seq.options <- seq.options

    x$paths <- data.frame(t(res$paths))
    rownames(x$paths) <- paste("X", 1:p, sep="")
    colnames(x$paths) <- NULL
    
  }
  
  # initial call to the function
  x$call <- match.call()
  
  class(x) <- "testHSIC"

  return(x)
  
}

### check.sensiHSIC ############################################################

check.sensiHSIC <- function(sensi, test.method){
  
  ### inputs ###
  
  # sensi:            1 x ...     list of class sensiHSIC (with all data)
  # test.method:      1 x 1       character (method used to test independence)
  
  ### output ###
  
  # --> If sensi and test.method are properly defined, the function returns NULL.
  # --> Otherwise, an error is returned indicating what happens.

  ### useful lists #############################################################
  
  L <- available.material()
  
  ### For sensi ################################################################
  
  msg <- "sensi is not an appropriate object of class sensiHSIC. \n"
  
#  if(!(class(sensi)=="sensiHSIC")){ # bug
   if(!(inherits(sensi,"sensiHSIC"))){
    stop("sensi must be a list of class sensiHSIC.", call.=FALSE)
  }
  
  ### For sensi$KX #############################################################
  
  # sensi must contain an element named KX
  if(is.null(sensi$KX)){
    stop(msg, "KX is missing in sensi. \n",
         "Please add an object named KX in sensi.", call.=FALSE)
  }else{
  
    msg.KX <- "sensi$KX must be a n-by-n-by-p array containing all p input Gram matrices."
    
    # sensi$KX must be an array
    if(!(is.array(sensi$KX))){
      stop(msg, msg.KX, call.=FALSE)
    }else{
      # sensi$KX must be a 3-dimensional array
      if(length(dim(sensi$KX))!=3){
        stop(msg, msg.KX, call.=FALSE)
      }else{
        n <- dim(sensi$KX)[1]
        p <- dim(sensi$KX)[3]
        # each slice of sensi$KX must be a square matrix
        if(!(ncol(sensi$KX)==n)) stop(msg, msg.KX, call.=FALSE)
      }
    }
    
  }
  
  ### For sensi$KY #############################################################
  
  # sensi must contain an element named KY
  if(is.null(sensi$KY)){
    stop(msg, "KY is missing in sensi. \n",
         "Please add an object named KY in sensi.", call.=FALSE)
  }else{
  
    msg.KY <- "sensi$KY must be a n-by-n array containing the output Gram matrix."
    
    # sensi$KY must be a matrix
    if(!(is.matrix(sensi$KY))){
      stop(msg, msg.KY, call.=FALSE)
    }else{
      # sensi$KX and sensi$KY must have equal dimensions
      if(nrow(sensi$KY)!=n | ncol(sensi$KY)!=n){
        stop(msg, "Dimensions of sensi$KX and sensi$KY do not match. \n",
             msg.KX, "\n", msg.KY, call.=FALSE)
      }
    }
    
  }

  ### For estimator.type #######################################################
  
  # sensi must contain an element named estimator.type
  if(is.null(sensi$estimator.type)){
    stop(msg, "estimator.type is missing in sensi. \n",
         "Please add an object named estimator.type in sensi.", call.=FALSE)
  }else{
  
    # sensi$estimator.type must be a string
    if(!is.character(sensi$estimator.type) | length(sensi$estimator.type)>1){
      stop(msg, "sensi$estimator.type must be a string.", call.=FALSE)
    }else{
      # sensi$estimator.type must be either "V-stat" or "U-stat"
      if(!(sensi$estimator.type %in% c("V-stat", "U-stat"))){
        stop(msg, "estimator.type must be: V-stat or U-stat.", call.=FALSE)
      }
    }
    
  }
  
  ### For cond #################################################################
  
  # sensi must contain an element named cond
  
  if(!is.null(sensi$cond)){
    
    # sensi$weights must then contain an element named weights
    if(is.null(sensi$weights)){
      stop(msg, "weights is missing in sensi. \n",
           "Please add an object named weights in sensi.", call.=FALSE)
    }else{

      msg.weights <- "sensi$weights must be a vector of length n."

      # weights must a vector
      if(!is.numeric(sensi$weights) && !is.null(dim(sensi$weights))){
        stop(msg, msg.weights, call.=FALSE)
      }else{
        # weights must be of length n
        if(length(sensi$weights)!=n){
          stop(msg, msg.weights, call.=FALSE)
        }
      }

    }
    
  }
  
  ### For test.method ##########################################################
  
  # in presence of U-statistics, test.method must belong to L$U.tests
  
  if(sensi$estimator.type=="U-stat"){
    if(!(test.method %in% L$U.tests)){
      stop("Invalid test method for U-statistics. \n",
           "In presence of U-statistics, possible test methods only include: ", 
           paste(L$U.tests, collapse=", "), ".", call.=FALSE)
    }
  }
  
  # in presence of C-HSIC indices, there are two constraints:
  # 1) cond$type must be "indicTh" (binary weights)
  # 2) test.method must belong to L$cond.tests
  
  if(!is.null(sensi$cond)){
    
    # cond$type must be "indicTh"
    if(!(sensi$cond$type=="indicTh")){
      stop("No available test method for C-HSIC indices ",
           "when sensi$cond$type is ", sensi$cond$type, ". \n",
           "sensi$cond$type=indicTh must be selected instead.", call.=FALSE)
    }
    
    # test.method must belong to L$cond.tests
    if(!(test.method %in% L$cond.tests)){
      stop("Invalid test method for C-HSIC indices. \n",
           "In presence of C-HSIC indices, possible test methods only include: ", 
           paste(L$cond.tests, collapse=", "), ".", call.=FALSE)
    }
    
  }
  
  ### For sensi$HSICXY and sensi$TO.num ########################################
  
  # two situations can be distinguished:
  # 1) the test is based on the numerators of the first-order HSIC-ANOVA indices
  # 2) __________________________________________ total-order __________________
  
  if(stringr::str_detect(test.method, "^Tot_")){
    
    # sensi must contain an element named TO.num
    if(is.null(sensi$TO.num)){
      stop(msg, "TO.num is missing in sensi. \n",
           "Please add an object named TO.num in sensi.", call.=FALSE)
    }else{
      
      # sensi$TO.num must be a data frame
      if(!is.data.frame(sensi$TO.num)){
        stop(msg, "sensi$TO.num must be a data frame.", call.=FALSE)
      }else{
        # sensi$TO.num must contain a column named "original"
        if(!("original" %in% colnames(sensi$TO.num))){
          stop(msg, "sensi$TO.num must have a column named original.", call.=FALSE)
        }else{
          # sensi$TO.num must have p rows
          if(nrow(sensi$TO.num)!=p){
            stop(msg, "Dimensions of sensi$KX and sensi$TO.num do not match. \n",
                 msg.KX, "\n", "sensi$TO.num must be a data frame with p rows.", call.=FALSE)
          }
        }
      }
      
    }
  
  }else{
    
    # sensi must contain an element named HSICXY
    if(is.null(sensi$HSICXY)){
      stop(msg, "HSICXY is missing in sensi. \n",
           "Please add an object named HSICXY in sensi.", call.=FALSE)
    }else{
      
      # sensi$HSICXY must be a data frame
      if(!is.data.frame(sensi$HSICXY)){
        stop(msg, "sensi$HSICXY must be a data frame.", call.=FALSE)
      }else{
        # sensi$HSICXY must contain a column named "original"
        if(!("original" %in% colnames(sensi$HSICXY))){
          stop(msg, "sensi$HSICXY must have a column named original.", call.=FALSE)
        }else{
          # sensi$HSICXY must have p rows
          if(nrow(sensi$HSICXY)!=p){
            stop(msg, "Dimensions of sensi$KX and sensi$HSICXY do not match. \n",
                 msg.KX, "\n", "sensi$HSICXY must be a data frame with p rows.", call.=FALSE)
          }
        }
      }
      
    }
 
  }
  
  return(NULL)
  
}

### asymptotic.test ############################################################

asymptotic.test <- function(HSIC.obs, KX, KY, estimator.type){
  
  ### inputs ###
  
  # HSIC.obs:         1 x p         numeric   (observed values of HSIC indices)
  # KX:               n x n x p     array     (input Gram matrices)
  # KY:               n x n         matrix    (output Gram matrix)
  # estimator.type:   1 x 1         character ("U-stat" or "V-stat")
  
  ### output ###
  
  # res:        1 x 2   list        (test results)
  
  # In particular, res contains the following objects:
  
  # pval:       1 x p   numeric     (estimated p-values)             
  # param:      p x 2   matrix      (parameters of the Gamma distributions)
  
  n <- dim(KX)[1]
  p <- dim(KX)[3]
  
  # pre-allocation
  pval <- rep(NA, p)
  param <- switch(estimator.type, 
                  "U-stat"={ matrix(NA, 3, p) },
                  "V-stat"={ matrix(NA, 2, p) }, NA)

  # output Gram matrix
  EY <- mean(diag(KY))            # mean of diagonal coefficients
  EYY <- mean.nondiag(KY)         # mean of non-diagonal coefficients
  BY <- double.centering(KY)      # the rows and columns are now centered
  
  for(i in 1:p){
    
    # input Gram matrix
    KXi <- KX[,,i]                  # i-th input Gram matrix
    EXi <- mean(diag(KXi))          # mean of diagonal coefficients
    EXiXi <- mean.nondiag(KXi)      # mean of non-diagonal coefficients
    BXi <- double.centering(KXi)    # double-centering 
    
    Bi <- (BXi*BY)^2
    
    # estimation of the mean and variance of the V-statistic
    HSIC.mean <- (EXi-EXiXi)*(EY-EYY)/n
    HSIC.var <- 2*(n-4)*(n-5)/(n*(n-1)*(n-2)*(n-3))*mean.nondiag(Bi)
    
    # method of moments
    alpha <- param[1,i] <- (HSIC.mean^2)/HSIC.var
    beta  <- param[2,i] <- HSIC.var/HSIC.mean
    if(estimator.type=="U-stat") delta <- param[3,i] <- -HSIC.mean
    
    # parametric estimation of the p-value
    pval[i] <- switch(estimator.type,
                      "U-stat"={ 1-pgamma(q=HSIC.obs[i]-delta, shape=alpha, scale=beta) },
                      "V-stat"={ 1-pgamma(q=HSIC.obs[i], shape=alpha, scale=beta) }, NA)
    
  }
  
  res <- list(pval=pval, param=param)
  
  return(res)
  
}

### HSIC.perm ##################################################################

HSIC.perm <- function(KX, KY, measure.def, estimator.type, B, A){
  
  ### inputs ###
  
  # KX:               n x n x p     array     (input Gram matrices)
  # KY:               n x n         matrix    (output Gram matrix)
  # measure.def:      1 x 1         character ("HSICXY" or "TO.num")
  # estimator.type:   1 x 1         character ("U-stat" or "V-stat")
  # B:                1 x 1         numeric   (number of permutations)
  # A:                1 x r         numeric   (subset of 1,...,p)
  
  ### output ###
  
  # Hperm is either a vector (if r=1) or a matrix (if r>1)
  
  # Hperm:            1 x B         numeric   (HSIC index after permutations)
  #                   B x r         matrix    (HSIC indices after permutations)
  
  n <- dim(KX)[1] # nb of samples
  p <- dim(KX)[3] # nb of input variables
  
  r <- length(A)
  Hperm <- matrix(NA, B, r)
  
  for(k in 1:r){
    
    i <- A[k]       # k-th element in the vector A
    KXi <- KX[,,i]  # i-th Gram matrix
    
    ### definition of the permutation scheme ###
      
    if(measure.def=="HSICXY"){

      perm.scheme <- function(j){ 
        shuf <- sample(seq(1,n), n)
        HSIC.fun(KXi, KY[shuf, shuf], estimator.type)
      }
      
    }else{
      
      kXi <- KXi-1      # the Gram matrix that needs to be permuted
      
      # product of all other Gram matrices
      LXi <- 1          
      for(s in setdiff(1:p, i)) LXi <- LXi*KX[,,s]
      # --> this matrix is invariant under permutations
      
      perm.scheme <- function(j){ 
        shuf <- sample(seq(1,n), n)
        HSIC.fun(kXi[shuf, shuf]*LXi, KY, estimator.type)
      }
      
    }

    ### repetition of B permutations ###
    
    Hperm[,k] <- sapply(1:B, perm.scheme)
    
  }
  
  if(r==1) Hperm <- as.numeric(Hperm)
  
  return(Hperm)
  
}

### permutation.test ###########################################################

permutation.test <- function(HSIC.obs, KX, KY,
                             measure.def, estimator.type, B){
  
  ### inputs ###
  
  # HSIC.obs:         1 x p         numeric   (observed values of HSIC indices)
  # KX:               n x n x p     array     (input Gram matrices)
  # KY:               n x n         matrix    (output Gram matrix)
  # measure.def:      1 x 1         character ("HSICXY" or "TO.num")
  # estimator.type:   1 x 1         character ("U-stat" or "V-stat")
  # B:                1 x 1         numeric   (number of permutations)
  
  ### output ###
  
  # res:        1 x 2   list        (test results)
  
  # In particular, res contains the following object:
  
  # pval:           1 x p   numeric     (estimated p-values)  
  # Hperm:          B x p   numeric     (HSIC indices after permutations)
  
  p <- dim(KX)[3]
  
  pval <- rep(NA, p)
  
  ### repetitions of B permutations ###
  
  Hperm <- HSIC.perm(KX, KY, measure.def, estimator.type, B, 1:p)
  
  ### non-parametric estimation of the p-values ###
  
  pval <- apply(mat.vec(Hperm, HSIC.obs, ">"), 2, mean)
  
  res <- list(pval=pval, Hperm=Hperm)
  
  return(res)
  
}

### seq.permutation.test #######################################################

seq.permutation.test <- function(HSIC.obs, KX, KY,
                                 measure.def, estimator.type, seq.options){
  
  ### inputs ###
  
  # HSIC.obs:         1 x p         numeric   (observed values of HSIC indices)
  # KX:               n x n x p     array     (input Gram matrices)
  # KY:               n x n         matrix    (output Gram matrix)
  # measure.def:      1 x 1         character ("HSICXY" or "TO.num")
  # estimator.type:   1 x 1         character ("U-stat" or "V-stat")
  # seq.options:      1 x ...       list      (options for the sequential test)
  
  ### output ###
  
  # res:        1 x 2     list      (test results)
  
  # In particular, res contains the following objects:
  
  # pval:       1 x p     numeric   (estimated p-values)
  # paths:    ... x p     matrix    (all estimated p-values)
  
  p <- dim(KX)[3]
  
  # hyperparameters to guide the sequential procedure
  criterion <- seq.options$criterion
  alpha <- seq.options$alpha
  Bstart <- seq.options$Bstart
  Bfinal <- seq.options$Bfinal
  Bbatch <- seq.options$Bbatch
  Bconv <- seq.options$Bconv
  graph <- seq.options$graph
  
  if(Bconv>Bstart){
    Bstart <- Bconv
  }

  ### goal-oriented test of independence ###
  
  if(criterion=="screening"){
    
    #############################
    ### FIRST CASE: SCREENING ###
    #############################
      
    pval <- rep(NA, p)              # p-values
    paths <- matrix(NA, Bfinal, p)  # estimation paths
    nb.perm <- rep(NA, p)           # nb of permutations
    
    for(i in 1:p){
      
      ### FIRST ITERATION ###
      
      B <- Bstart # counter for permutations
      
      # observed value of the HSIC index
      Hi <- HSIC.obs[i]
      # new values of the HSIC index after Bstart permutations 
      seq.HSIC <- HSIC.perm(KX, KY, measure.def, estimator.type, B, i)
      # estimation of the p-value with increasing sample sizes
      seq.pval <- cumsum(seq.HSIC>Hi)/seq(2, B+1)
      # screening (0 if pval>alpha, 1 otherwise)
      seq.bin <- 1*(seq.pval<=alpha)
      # more permutations?
      more.perm <- !identical( seq.bin[(B-Bconv+1):B], 
                               rep(seq.bin[B], Bconv) )
      
      ### OTHER ITERATIONS ###
      
      while(more.perm && B<Bfinal){
        
        B <- B+Bbatch
        
        # new values of the HSIC index (coming from Bbatch additional permutations)
        seq.HSIC <- c(seq.HSIC,
                      HSIC.perm(KX, KY, measure.def, estimator.type, Bbatch, i))
        # all p-values are recomputed
        seq.pval <- cumsum(seq.HSIC>Hi)/seq(2, B+1)
        # all binary decisions are recomputed
        seq.bin <- 1*(seq.pval<=alpha)       
        # more permutations?
        more.perm <- !identical( seq.bin[(B-Bconv+1):B], 
                                 rep(seq.bin[B], Bconv) )
        
      }
      
      # storage
      pval[i] <- seq.pval[B]
      paths[1:B,i] <- seq.pval
      nb.perm[i] <- B
      
    }
    
    Bm <- max(nb.perm)    # largest nb of permutations
    paths <- paths[1:Bm,] # paths are shortened
    
  }else{
    
    ####################################
    ### SECOND CASE: RANKING or BOTH ###
    ####################################
    
    ### FIRST ITERATION ###
    
    B <- Bstart # counter for permutations
    
    # new values of the HSIC indices after Bstart permutations 
    seq.HSIC <- HSIC.perm(KX, KY, measure.def, estimator.type, B, 1:p)
    # estimation of the p-value with increasing sample sizes
    seq.pval <- ( apply(mat.vec(seq.HSIC, HSIC.obs, ">"), 2, cumsum) /
                    matrix(seq(2, B+1), B, p) )
    # rankings of p-values with increasing sample sizes
    seq.rk <- t(apply(seq.pval, 1, rank, ties="first"))
    # more permutations?
    more.perm <- !identical( seq.rk[(B-Bconv+1):B,], 
                             matrix(seq.rk[B,], Bconv, p, byrow=TRUE) )
    
    # additional check if screening="both"
    if(!more.perm && criterion=="both"){
      # screening (0 if pval>alpha, 1 otherwise)
      seq.bin <- 1*(seq.pval<=alpha)
      # more permutations?
      more.perm <- !identical( seq.bin[(B-Bconv+1):B,], 
                               matrix(seq.bin[B,], Bconv, p, byrow=TRUE) )
    }
    
    ### OTHER ITERATIONS ###
    
    while(more.perm && B<Bfinal){
      
      B <- B+Bbatch
      
      # new values of the HSIC indices (coming from Bbatch additional permutations)
      seq.HSIC <- rbind(seq.HSIC,
                        HSIC.perm(KX, KY, measure.def, estimator.type, Bbatch, 1:p))
      # all p-values are recomputed
      seq.pval <- ( apply(mat.vec(seq.HSIC, HSIC.obs, ">"), 2, cumsum) /
                      matrix(seq(2, B+1), B, p) )
      # all rankings are recomputed
      seq.rk <- t(apply(seq.pval, 1, rank, ties="first"))
      # more permutations?
      more.perm <- !identical( seq.rk[(B-Bconv+1):B,], 
                               matrix(seq.rk[B,], Bconv, p, byrow=TRUE) )
      
      # additional check if screening="both"
      if(!more.perm && criterion=="both"){
        seq.bin <- 1*(seq.pval<=alpha)       
        # more permutations?
        more.perm <- !identical( seq.bin[(B-Bconv+1):B,], 
                                 matrix(seq.bin[B,], Bconv, p, byrow=TRUE) )
      }
      
    }
    
    # storage
    pval <- seq.pval[B,]
    paths <- seq.pval
    
  }
  
  ### warning if paths are not stabilized ###
  
  if(criterion=="screening"){
    
    if(any(nb.perm>Bfinal)){
      ind <- which(nb.perm>Bfinal)
      warning("Even after ", Bfinal, " iterations, ",
              "the algorithm has not converged yet for the following variables: ",
              paste(ind, collapse=", "), ".", call.=FALSE)
    }
    
  }else{
    
    if(B>Bfinal){
      warning("Even after ", Bfinal, " iterations, ",
              "the algorithm has not converged yet.", call.=FALSE)
    }
    
  }
  
  ### plot of estimation paths ###
  
  if(graph){
    # estimation paths for all p-values
    matplot(paths, type='l', ylim=c(0,1),
            xlab="number of permutations",
            ylab="estimated p-values",
            main="Sequential estimation of the p-values until convergence")
    # red line for screening
    if(!(criterion=="ranking")) abline(h=alpha, lwd=2, col="red")
  }
  
  res <- list(pval=pval, paths=paths)
  
  return(res)
  
}

### gamma.test #################################################################

gamma.test <- function(HSIC.obs, KX, KY, measure.def){

  ### inputs ###

  # HSIC.obs:       1 x p         numeric   (observed values of HSIC indices)
  # KX:             n x n x p     array     (input Gram matrices)
  # KY:             n x n         matrix    (output Gram matrix)
  # measure.def:    1 x 1         character ("HSICXY" or "TO.num")

  ### output ###

  # res:        1 x 2   list        (test results)

  # In particular, res contains the following objects:

  # pval:       1 x p   numeric     (estimated p-values)
  # param:      2 x p   matrix      (parameters of the Gamma distributions)

  # pre-allocation
  
  n <- dim(KX)[1]
  p <- dim(KX)[3]
  pval <- rep(NA, p)
  param <- matrix(NA, 2, p)
  
  # computations of the matrices Ai and Wi for all i in {1,...,p}
  
  mat <- compute.mat.TrAW(KX, KY, measure.def)
  
  # computations of the first two moments of Tr(Ai Wi) for all i in {1,...,p}
  
  mom <- compute.mom.TrAW(mat$A, mat$W, measure.def)
  
  # parametric estimation of all p-values
  
  for(i in 1:p){
    
    esp <- mom[1,i]
    var <- mom[2,i]

    # method of moments for Tr(Ai W)
    shape.TrAW <- (esp^2)/var
    scale.TrAW <- var/esp
    
    # parameters for HSIC(Xi,Y) = Tr(Ai W)/(n^2)
    alpha <- param[1,i] <- shape.TrAW
    beta  <- param[2,i] <- scale.TrAW/(n^2)
    
    # parametric estimation of the p-value
    pval[i] <- 1-pgamma(q=HSIC.obs[i], shape=alpha, scale=beta)
    
  }
  
  res <- list(pval=pval, param=param)

  return(res)

}

### compute.mat.TrAW ###########################################################

compute.mat.TrAW <- function(KX, KY, measure.def){
  
  ### inputs ###
  
  # KX:             n x n x p     array     (input Gram matrices)
  # KY:             n x n         matrix    (output Gram matrix)
  # measure.def:    1 x 1         character ("HSICXY" or "TO.num")
  
  ### output ###
  
  # mat:      1 x 2         list      (all matrices)
  
  # In particular, mat contains the following objects:
  
  # A:        n x n x p     array     (all matrices Ai)
  
  # W:        n x n         matrix    (matrix W        if measure.def="HSICXY")
  #           n x n x p     array     (all matrices Wi if measure.def="TO.num")
  
  n <- dim(KX)[1]
  p <- dim(KX)[3]
  
  if(measure.def=="HSICXY"){
    
    A <- array(NA, c(n, n, p))
    for(i in 1:p) A[,,i] <- double.centering(KX[,,i])
    
    W <- double.centering(KY)

  }else{
    
    A <- array(NA, c(n, n, p))
    LY <- double.centering(KY)
    for(i in 1:p){
      Ai <- LY
      for(j in setdiff(1:p, i)) Ai <- Ai*KX[,,j]
      A[,,i] <- Ai
    }
    
    W <- KX-1
    
  }
  
  mat <- list(A=A, W=W)
  return(mat)
  
}

### compute.mom.TrAW ###########################################################

compute.mom.TrAW <- function(A, W, measure.def){
  
  ### inputs ###
  
  # A:              n x n x p     array     (all matrices Ai)
  
  # W:              n x n         matrix    (matrix W        if measure.def="HSICXY")
  #                 n x n x p     array     (all matrices Wi if measure.def="TO.num")
  
  # measure.def:    1 x 1         character ("HSICXY" or "TO.num")

  ### output ###
  
  # mom:    2 x p   numeric   (expectations and variances of Tr(Ai Wi) for all i in {1,...,p})  
  
  n <- dim(A)[1]
  p <- dim(A)[3]
  
  mom <- matrix(NA, 2, p)
  
  if(measure.def=="HSICXY"){
    
    # denominators
    denom1 <- ((n-1)^2)*(n+1)*(n-2)
    denom2 <- (n+1)*n*(n-1)*(n-2)*(n-3)
    
    # matrix operations
    tr.W <- sum(diag(W))      # complexity:          1 x n 
    tr.W2 <- sum(diag(W)^2)   # complexity:          2 x n
    sum.W2 <- sum(W^2)        # complexity: 2 x n2
                              # __________________________
                              # complexity: 2 x n2 + 3 x n
    
    # terms to be used in the final formulas
    O1.W <- (n-1)*sum.W2-tr.W^2
    O2.W <- n*(n+1)*tr.W2-(n-1)*(tr.W^2+2*sum.W2)
    
    for(i in 1:p){
      
      Ai <- A[,,i]
      
      # matrix operations
      tr.Ai <- sum(diag(Ai))
      tr.Ai2 <- sum(diag(Ai)^2)
      sum.Ai2 <- sum(Ai^2)
      
      # terms to be used in the final formulas
      O1.Ai <- (n-1)*sum.Ai2-tr.Ai^2
      O2.Ai <- n*(n+1)*tr.Ai2-(n-1)*(tr.Ai^2+2*sum.Ai2)
      
      # --> final formulas
      
      mom[1,i] <- tr.Ai*tr.W/(n-1)
      mom[2,i] <- 2*O1.Ai*O1.W/denom1 + O2.Ai*O2.W/denom2

    }

  }else{
    
    # denominators
    denom1 <- n
    denom2 <- n*(n-1)
    denom3 <- n*(n-1)*(n-2)
    denom4 <- n*(n-1)*(n-2)*(n-3)
    
    for(i in 1:p){

      # --> for the matrix Ai
      
      Ai <- A[,,i]
      tAi <- Ai
      diag(tAi) <- 0
      
      # matrix operations
      tr.Ai <- sum(diag(Ai))                      # complexity:          1 x n
      sum.tAi <- sum(tAi)                         # complexity: 1 x n2
      tr.Ai2 <- sum(diag(Ai)^2)                   # complexity:          2 x n
      sum.Ai2 <- sum(Ai*Ai)                       # complexity: 2 x n2
      sdc.Ai <- sum(diag(Ai)*colSums(Ai))         # complexity: 1 x n2 + 2 x n
      sum.tAi2 <- sum.Ai2-tr.Ai2          
      sct.Ai <- sum(colSums(tAi)%*%tAi)-sum.tAi2  # complexity: 3 x n2 + 1 x n
                                                  # __________________________
                                                  # complexity: 7 x n2 + 6 x n 
      
      # terms to be used in the final formulas
      O1.Ai <- tr.Ai2
      O21.Ai <- (tr.Ai)^2-tr.Ai2
      O22.Ai <- sum.tAi2
      O23.Ai <- sdc.Ai-tr.Ai2
      O31.Ai <- sct.Ai
      O32.Ai <- tr.Ai*sum.tAi-2*(sdc.Ai-tr.Ai2)
      O4.Ai <- (sum.tAi)^2-2*sum.tAi2-4*sct.Ai
      
      # --> for the matrix Wi
      
      Wi <- W[,,i]
      tWi <- Wi
      diag(tWi) <- 0
      
      # matrix operations
      tr.Wi <- sum(diag(Wi))
      sum.tWi <- sum(tWi)
      tr.Wi2 <- sum(diag(Wi)^2)
      sum.Wi2 <- sum(Wi*Wi)
      sdc.Wi <- sum(diag(Wi)*colSums(Wi))
      sum.tWi2 <- sum.Wi2-tr.Wi2
      sct.Wi <- sum(colSums(tWi)%*%tWi)-sum.tWi2
      
      # terms to be used in the final formulas
      O1.Wi <- tr.Wi2
      O21.Wi <- (tr.Wi)^2-tr.Wi2
      O22.Wi <- sum.tWi2
      O23.Wi <- sdc.Wi-tr.Wi2
      O31.Wi <- sct.Wi
      O32.Wi <- tr.Wi*sum.tWi-2*(sdc.Wi-tr.Wi2)
      O4.Wi <- (sum.tWi)^2-2*sum.tWi2-4*sct.Wi
      
      # crossed expressions
      T1 <- ( O1.Ai*O1.Wi )/denom1
      T2 <- ( O21.Ai*O21.Wi + 2*O22.Ai*O22.Wi + 4*O23.Ai*O23.Wi )/denom2
      T3 <- ( 4*O31.Ai*O31.Wi + 2*O32.Ai*O32.Wi )/denom3
      T4 <- ( O4.Ai*O4.Wi )/denom4
      
      # --> final formulas
      
      E1 <- tr.Ai*tr.Wi/n+sum.tAi*sum.tWi/(n*(n-1))
      E2 <- T1+T2+T3+T4
      mom[, i] <- c(E1, E2-E1^2)
 
    }

  }
  
  return(mom)
  
}

### mean.nondiag ###############################################################

mean.nondiag <- function(A){
  
  ### input ###
  
  # A:    n x n     numeric (matrix)
  
  ### output ###
  
  # s:    1 x 1     numeric (mean of non-diagonal coefficients)
  
  n <- nrow(A)
  s <- sum(A-diag(diag(A)))/(n*(n-1))
  
  return(s)
  
}

### double.centering ###########################################################

double.centering <- function(A){
  
  ### input ###
  
  # A:        n x m         matrix    (initial given matrix)
  
  ### output ###
  
  # B:        n x m         matrix    (double-centered matrix)
  
  n <- nrow(A)
  m <- ncol(A)
  
  S.mat <- sum(A)/(n*m)
  S.rows <- rowSums(A)/m
  S.cols <- colSums(A)/n
  
  B <- ( A - matrix(rep(S.rows, m), n, m)
         - matrix(rep(S.cols, n), n, m, byrow=TRUE) + S.mat )
  
  return(B)
  
}

### mat.vec ####################################################################

mat.vec <- function(M1, vec, sgn){
  
  ### input ###
  
  # M1:     n x p     matrix      (samples from a random variable in R^p)  
  # vec:    1 x p     numeric     (reference vector in R^p)
  # sgn:    1 x 1     character   (boolean operator provided as a string)
  
  ### output ###
  
  # M2:     n x p     matrix      (each row in M1 is compared to vec with sgn)
  
  n <- nrow(M1)
  
  vec <- t(as.matrix(vec))
  vec.cloned <- matrix(1, n, 1)%*%vec
  M2 <- switch(sgn, 
               "==" ={ 1*(M1-vec.cloned==0) },
               ">"  ={ 1*(M1-vec.cloned>0) },
               "<"  ={ 1*(M1-vec.cloned<0) },
               ">=" ={ 1*(M1-vec.cloned>=0) },
               "<=" ={ 1*(M1-vec.cloned<=0) },
               NA)
  
  return(M2)
  
}

### print.testHSIC #############################################################

print.testHSIC <- function(x, ...){
  
  ### input ###
  
  # x:      1 x 1     list of class testHSIC  
  
  ##########
  
  # initial call
  cat("\n", "Call:", "\n", deparse(x$call), "\n\n", sep="")
  
  if(!is.null(x$pval)){
    
    # p-values
    cat("P-values:", "\n", sep="")
    print(x$pval)
    cat("\n")
    
    # parametric distributions
    if(x$prop[2]=="parametric"){
      cat("The test statistics were assumed to follow", x$family, "distributions:", "\n")
      print(x$param)
      cat("\n")
    }
    
  }else{
    cat("(empty)", "\n\n")
  }

}

### plot.testHSIC ##############################################################

plot.testHSIC <- function(x, ylim=c(0, 1), err=NA, ...){
  
  ### inputs ###
  
  # x:      1 x 1     list of class sensiHSIC 
  # ylim:   1 x 2     (numeric)                 bounds on the y-axis
  # err:    1 x 1     (numeric)                 level of type I error
  
  ##########
  
  if(is.null(x$pval)){
    
    warning("No result is stored in x. \n", call.=TRUE)
    
  }else{
    
    ###################
    ### Preparation ###
    ###################
    
    yshift <- 0.2
    
    # minimum and maximum values on the y-axis
    ym <- ylim[1]-yshift
    yM <- ylim[2]+2*yshift
    
    # nb of input variables
    p <- nrow(x$pval)
    
    # type I error
    if(is.na(err)){
      if(!is.null(x$seq.options)){
        err <- x$seq.options$alpha
      }else{
        err <- 0.05
        warning("By default, err=0.05 has been selected.", call.=FALSE)
      }
    }
    
    # preliminary work for the legend
    
    cap.H0 <- "Rejected variables"
    cap.H1 <- "Selected variables"
    
    col.H0 <- "red"
    col.H1 <- "green"
    
    all.captions <- c(cap.H0, cap.H1)
    all.colors <- c(col.H0, col.H1)
    
    #########################################
    ### Rejected variables (decision: H0) ###
    #########################################
    
    if(sum(x$pval>=err)>=1){
      var.H0 <- which(x$pval>=err)
      nodeplot(x$pval[var.H0, "original", drop=FALSE], at=var.H0, 
               xlim=c(1,p), ylim=c(ym, yM), pch=21, bg=col.H0)
    }
    
    #########################################
    ### Selected variables (decision: H1) ###
    #########################################
    
    if(sum(x$pval<err)>=1){
      var.H1 <- which(x$pval<err)
      nodeplot(x$pval[var.H1, "original", drop=FALSE], at=var.H1, 
               xlim=c(1,p), ylim=c(ym, yM), pch=21, bg=col.H1,
               add=TRUE)
    }
    
    ##############
    ### Legend ###
    ##############
    
    abline(h=err, col="red", lwd=2)
    abline(h=0, col="black", lwd=2, lty=2)
    abline(h=1, col="black", lwd=2, lty=2)
    
    legend("top", legend=all.captions, col="black", pch=21, pt.bg=all.colors)
    
  }
  
}