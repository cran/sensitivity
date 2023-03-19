

sensiHSIC <- function(model = NULL, X, target = NULL, cond = NULL, 
                      kernelX = "rbf", paramX = NA,
                      kernelY = "rbf", paramY = NA,
                      estimator.type = "V-stat",
                      nboot = 0, conf = 0.95,
                      anova = list(obj = "no", is.uniform = TRUE),
                      sensi = NULL, 
                      save.GM = list(KX = TRUE, KY = TRUE), ...){
  
  ### useful lists #############################################################
  
  L <- available.material()

  ### For model ################################################################
  
  # model must be a function
  
  if(!is.null(model)){
    if(!is.function(model)){
      stop("model must be a function.", call.=FALSE)
    }
  }
  
  ### For X ####################################################################
  
  if(missing(X)){
    stop("You must specificy what X is.", call.=FALSE)
  }
  
  # X must be a matrix or a data frame
  
  if(is.data.frame(X)){
    X <- as.matrix(unname(X))
  }else{
    if(!is.matrix(X)){
      stop("X must be a matrix or a data frame.", call.=FALSE)
    }
  }
  
  ### For model and X ##########################################################
  
  # if model is a function, it can be applied on the first row to compute the first output
  
  if(is.function(model)){
    tryCatch(model(X[1,,drop=FALSE]),
             error=function(e){
               stop("model is not designed to compute outputs from X. \n",
                    "Remember that model is applied to each row vector in X.", call.=FALSE)
               return(NA) },
             warning=function(w){
               return(NULL) }
             )
  }
  
  ### For target ###############################################################
  
  if(!is.null(target)){
    
    ### minimum requirements ###
    
    # 1) target must be a list of options
    # 2) target must contain at least one element named "c" (for the threshold)
  
    if(!is.list(target)){
      stop("target must be a list of options.", call.=FALSE)
    }else{
      if(is.null(target$c)){
        stop("Threshold not defined. You must specify target$c.", call.=FALSE)
      }
    }
  
    ### default values for other options ###
    
    if(is.null(target$type)) target$type <- "indicTh"
    if(is.null(target$upper)) target$upper <- TRUE
    if(is.null(target$param)) target$param <- 1
    
    ### removal of undesired options in target ###
    
    expected.target <- c("c", "type", "upper", "param")
    found.target <- names(target)
    if(!(all(found.target %in% expected.target))){
      undesired.target <- setdiff(found.target, expected.target)
      target[undesired.target] <- NULL
      warning("The list target can only contain the following options: ",
              paste0(expected.target, collapse=', '), ". \n", 
              "Other options have been removed.", call.=FALSE)
    }
    
    ### check of remaining options in target ###
    
    # target$c must a real number
    
    if(!is.numeric(target$c) | length(target$c)>1){
      stop("target$c must be a numeric of length 1.", call.=FALSE)
    }
    
    # target$upper must be a boolean 
    
    if(length(target$upper)>1){
      stop("target$upper must be a boolean.", call.=FALSE)
    }else{
      # must be either TRUE or FALSE
      if(!(isTRUE(target$upper) | isFALSE(target$upper))){
        stop("target$upper must be a boolean.", call.=FALSE)
      }
    }
    
    # target$type must be a string in L$weight.functions
    
    if(!is.character(target$type) | length(target$type)>1){
      stop("target$type must be a string.", call.=FALSE)
    }else{
      if(!(target$type %in% L$weight.functions)){
        target$type <- "exp1side"
        warning("target$type must be: indicTh, zeroTh, logistic or exp1side. \n",
                "By default, target$type=exp1side has been selected.", call.=FALSE)
      }
    }
    
    # target$param must be a positive real number (or NA)
    
    # if target$param=NA, default assignment
    if(length(target$param)==1){
      if(is.na(target$param)){
        target$param <- 1
      }
    }
    
    if(!is.numeric(target$param) | length(target$param)>1){
      stop("target$param must be a numeric of length 1.", call.=FALSE)
    }else{
      if(target$param<=0){
        target$param <- 1
        warning("target$param must be a positive real number. \n",
                "By default, target$param=1 has been selected.", call.=FALSE)
      }
    }
    
    # --> see tell.sensiHSIC for the checks based on output samples
    
  }
  
  ### For cond #################################################################
  
  if(!is.null(cond)){
    
    ### minimum requirements ###
    
    # 1) cond must be a list of options
    # 2) cond must contain at least one element named "c" (for the threshold)
    
    if(!is.list(cond)){
      stop("cond must be a list of options.", call.=FALSE)
    }else{
      if(is.null(cond$c)){
        stop("Threshold not defined. You must specify cond$c.", call.=FALSE)
      }
    }
    
    ### default values for other options ###
    
    if(is.null(cond$type)) cond$type <- "exp1side"
    if(is.null(cond$upper)) cond$upper <- TRUE
    if(is.null(cond$param)) cond$param <- 1
    
    ### removal of undesired options in cond ###
    
    expected.cond <- c("c", "type", "upper", "param")
    found.cond <- names(cond)
    if(!(all(found.cond %in% expected.cond))){
      undesired.cond <- setdiff(found.cond, expected.cond)
      cond[undesired.cond] <- NULL
      warning("The list cond can only contain the following options: ",
              paste0(expected.cond, collapse=', '), ". \n", 
              "Other options have been removed.", call.=FALSE)
    }
    
    ### check of remaining options in cond ###
    
    # cond$c must be a real number
    
    if(!is.numeric(cond$c) | length(cond$c)>1){
      stop("cond$c must be a numeric of length 1.", call.=FALSE)
    }
    
    # cond$upper must be a boolean 
    
    if(length(cond$upper)>1){
      stop("cond$upper must be a boolean.", call.=FALSE)
    }else{
      # must be either TRUE or FALSE
      if(!(isTRUE(cond$upper) | isFALSE(cond$upper))){
        stop("cond$upper must be a boolean.", call.=FALSE)
      }
    }
    
    # cond$type must be a string in L$weight.functions
    
    if(!is.character(cond$type) | length(cond$type)>1){
      stop("cond$type must be a string.", call.=FALSE)
    }else{
      if(!(cond$type %in% L$weight.functions)){
        cond$type <- "exp1side"
        warning("cond$type must be: indicTh, zeroTh, logistic or exp1side. \n",
                "By default, cond$type=exp1side has been selected.", call.=FALSE)
      }
    }
    
    # cond$param must be a positive real number (or NA)
    
    # if cond$param=NA, default assignment
    if(length(cond$param)==1){
      if(is.na(cond$param)){
        cond$param <- 1
      }
    }
    
    if(!is.numeric(cond$param) | length(cond$param)>1){
      stop("cond$param must be a numeric of length 1.", call.=FALSE)
    }else{
      if(cond$param<=0){
        cond$param <- 1
        warning("cond$param must be a positive real number. \n",
                "By default, cond$param=1 has been selected.", call.=FALSE)
      }
    }
    
    # --> see tell.sensiHSIC for the checks based on output samples
    
  }
  
  ### For target and cond ######################################################
  
  # if target and cond are enabled simultaneously, cond is ignored
  
  if(!is.null(target) && !is.null(cond)){
    cond <- NULL
    warning("target and cond cannot be enabled simultaneously. \n",
            "cond has been ignored.", call.=FALSE)
  }
  
  ### For kernelX ##############################################################
  
  p <- ncol(X)
  nkX <- length(kernelX)
  
  ### minimum requirements ###
  
  # 1) kernelX must be a string or a vector of strings
  # 2) kernelX must be of length 1 or p
  # 3) The string(s) in kernelX must all belong to L$kernels
  
  if(!is.character(kernelX)){
    stop("kernelX must be a string or a vector of strings.", call.=FALSE)
  }
  
  if(!(nkX==1 | nkX==p)){
    stop("kernelX must be of length 1 or p (number of input variables).", call.=FALSE)
  }
  
  if(!all(kernelX %in% L$kernels)){
    stop("Invalid kernel names in kernelX. \n", 
         "Available kernels include: ",
         paste0(L$kernels, collapse=', '), ".", call.=FALSE)
  }

  ### For paramX ###############################################################
  
  ### minimum requirements ###
  
  # 1) paramX must be of length 1 or p
  # 2) paramX must be a scalar value or a vector of values
  # 3) elements in paramX must be NA or numeric
  
  npX <- length(paramX)
  
  if(!(npX==1 | npX==p)){
    stop("paramX must be of length 1 or p (number of input variables).", call.=FALSE)
  }
  
  msg <- "paramX must be a scalar value or a vector of values."
  
  if(is.numeric(paramX)){
    
    if(!is.null(dim(paramX))) stop(msg, call.=FALSE)
    
  }else if(is.logical(paramX)){
    
    if(!all(is.na(paramX))) stop(msg, call.=FALSE)
    
  }else{
    
    stop(msg, call.=FALSE)
    
  }

  ### For kernelY ##############################################################
  
  ### minimum requirements ###
  
  # 1) kernelY can be either:
  # --> a) a string or a vector of strings 
  # --> b) a list of options
  
  # If a): 
  # --> the string(s) in kernelY must all belong to L$kernels
  
  # If b):
  # --> kernelY must contain at least one element named "method"
  # --> kernelY$method must be: PCA, GAK or DTW
  
  if(!is.character(kernelY) && !is.list(kernelY)){
    stop("kernelY must be a string, a vector of strings or a list of options.", call.=FALSE)
  }else{
    if(is.character(kernelY)){
    # kernelY is a string or a vector of strings
      if(!all(kernelY %in% L$kernels)){
        stop("Invalid kernel names in kernelY. \n", 
             "Available kernels include: ",
             paste0(L$kernels, collapse=', '), ".", call.=FALSE)
      }
    }else{
    # kernelY is a list of options
      # must contain an object named "method"
      if(is.null(kernelY$method)){
        stop("Method not defined. You must specify kernelY$method.", call.=FALSE)
      }else{
        # kernelY must be a string
        if(!is.character(kernelY$method) | length(kernelY$method)>1){
          stop("kernelY$method must be a string.", call.=FALSE)
        }else{
          # must be: PCA, GAK or DTW
          list.kY <- c("PCA", "GAK", "DTW")
          if(!(kernelY$method %in% list.kY)){
            stop("kernelY$method must be: ",
                 paste0(list.kY, collapse=', '), ".", call.=FALSE)
          }
        }
      }
    }
  }

  ### If b): default values and other checks ###

  if(is.list(kernelY)){

    if(kernelY$method=="PCA"){
      
      ### default values for missing options ###
      
      if(is.null(kernelY$data.centering)) kernelY$data.centering <- TRUE
      if(is.null(kernelY$data.scaling)) kernelY$data.scaling <- TRUE
      if(is.null(kernelY$fam)) kernelY$fam <- "rbf"
      if(is.null(kernelY$combi)) kernelY$combi <- "sum"
      if(is.null(kernelY$position)) kernelY$position <- "intern"
      
      ### default values for expl.var and PC ###
      
      # 1) default values 
      if(is.null(kernelY$expl.var)){
        if(is.null(kernelY$PC)){ 
          # expl.var is NULL + PC is NULL
          kernelY$expl.var <- 0.95
          kernelY$PC <- NA
        }else{
          # expl.var is NULL + PC has a value
          kernelY$expl.var <- NA
        }
      }else{
        if(is.null(kernelY$PC)){
          # expl.var has a value + PC is NULL
          kernelY$PC <- NA
        }
      }
      
      # 2) lack of information / excess of information
      expl.var.NA <- isTRUE(is.na(kernelY$expl.var))
      PC.NA <- isTRUE(is.na(kernelY$PC))
      # three situations:
      # a) if the two values are NA --> set kernelY$expl.var=0.95
      # b) if none of the two values is NA --> set kernelY$PC=NA
      # c) if only one of the two values is NA --> do nothing
      if(expl.var.NA && PC.NA){
        kernelY$expl.var <- 0.95
        warning("kernelY$expl.var and kernelY$PC cannot be both equal to NA. \n",
                "By default, kernelY$expl.var=0.95 has been selected.", call.=FALSE)
      }else{
        if(!expl.var.NA && !PC.NA){
          kernelY$PC <- NA
          warning("kernelY$expl.var and kernelY$PC were enabled simultaneously. \n",
                  "kernelY$PC has been ignored.", call.=FALSE)
        }
      }
      
      ### removal of undesired options in kernelY ###
      
      expected.kY <- c("method", 
                       "data.centering", "data.scaling",
                       "fam", "expl.var", "PC", "combi", "position")
      found.kY <- names(kernelY)
      if(!(all(found.kY %in% expected.kY))){
        undesired.kY <- setdiff(found.kY, expected.kY)
        kernelY[undesired.kY] <- NULL
        warning("The list kernelY can only contain the following options: ",
                paste0(expected.kY, collapse=', '), ". \n",
                "Other options have been removed.", call.=FALSE)
      }
      
      ### check of remaining options in kernelY ###
      
      # kernelY$data.centering must be a boolean
      
      if(length(kernelY$data.centering)>1){
        stop("kernelY$data.centering must be a boolean.", call.=FALSE)
      }else{
        if(!(isTRUE(kernelY$data.centering) | isFALSE(kernelY$data.centering))){
          stop("kernelY$data.centering must be a boolean.", call.=FALSE)
        }
      }
      
      # kernelY$data.scaling must be a boolean
      
      if(length(kernelY$data.scaling)>1){
        stop("kernelY$data.scaling must be a boolean.", call.=FALSE)
      }else{
        if(!(isTRUE(kernelY$data.scaling) | isFALSE(kernelY$data.scaling))){
          stop("kernelY$data.scaling must be a boolean.", call.=FALSE)
        }
      }

      # kernelY$fam must be a string.
      # the string in kernelY$fam must belong to L$pca.kernels
      
      if(!is.character(kernelY$fam) | length(kernelY$fam)>1){
        stop("kernelY$fam must be a string.", call.=FALSE)
      }else{
        if(!(kernelY$fam %in% L$pca.kernels)){
          stop("Invalid kernel name in kernelY$fam. \n", 
               "Available kernels include: ",
               paste0(L$pca.kernels, collapse=', '), ".", call.=FALSE)
        }
      }
      
      # kernelY$expl.var must be real number between 0 and 1 (or NA)
      
      if(length(kernelY$expl.var)>1){
        stop("kernelY$expl.var must be a real number between 0 and 1.", call.=FALSE)
      }else{
        if(!is.numeric(kernelY$expl.var) && !is.na(kernelY$expl.var)){
          stop("kernelY$expl.var must be a real number between 0 and 1.", call.=FALSE)
        }else{
          if(!is.na(kernelY$expl.var)){
            if(kernelY$expl.var<0 | kernelY$expl.var>1){
              kernelY$expl.var <- 0.95
              warning("kernelY$expl.var was not between 0 and 1. \n",
                      "0.95 has been taken instead.", call.=FALSE)
            }
          }
        }
      }

      # kernelY$PC must be a positive integer (or NA)
      
      if(length(kernelY$PC)>1){
        stop("kernelY$PC must be a positive integer.", call.=FALSE)
      }else{
        if(!is.numeric(kernelY$PC) && !is.na(kernelY$PC)){
          stop("kernelY$PC must be a positive integer.", call.=FALSE)
        }else{
          if(!is.na(kernelY$PC)){
            if(kernelY$PC<=0){
              stop("kernelY$PC must be a positive integer.", call.=FALSE)
            }else{
              if(!(kernelY$PC%%1==0)){
                kernelY$PC <- floor(kernelY$PC)
                warning("kernelY$PC was not an integer. \n",
                        "floor(kernelY$PC) has been taken instead.", call.=FALSE)
              }
            }
          }
        }
      }
 
      # kernelY$combi must be either "sum" or "prod"
      
      if(!is.character(kernelY$combi) | length(kernelY$combi)>1){
        stop("kernelY$combi must be a string.", call.=FALSE)
      }else{
        if(!(kernelY$combi %in% c("sum", "prod"))){
          kernelY$combi <- "sum"
          warning("kernelY$combi must be: sum or prod. \n",
                  "By default, kernelY$combi=sum has been selected.", call.=FALSE)
        }
      }
      
      # kernelY$position must be either "no", "intern" or "extern"
      
      if(!is.character(kernelY$position) | length(kernelY$position)>1){
        stop("kernelY$position must be a string.", call.=FALSE)
      }else{
        if(!(kernelY$position %in% c("nowhere", "intern", "extern"))){
          kernelY$position <- "intern"
          warning("kernelY$position must be: nowhere, intern or extern. \n",
                  "By default, kernelY$position=intern has been selected.", call.=FALSE)
        }
      }

      # if kernelY$combi="prod", kernelY$position must be "intern" 
      
      if(kernelY$combi=="prod" && kernelY$position=="extern"){
        kernelY$position <- "intern"
        warning("If kernelY$combi=prod, kernelY$position cannot be extern. \n",
                "By default, kernelY$position=intern has been selected.", call.=FALSE)
      }

    }
    
    if(kernelY$method=="GAK"){
      
      ### removal of undesired options in kernelY ###
      
      expected.kY <- "method"
      found.kY <- names(kernelY)
      if(!(all(found.kY %in% expected.kY))){
        undesired.kY <- setdiff(found.kY, expected.kY)
        kernelY[undesired.kY] <- NULL
        warning("The list kernelY can only contain the following option: method. \n",
                "Other options have been removed.", call.=FALSE)
      }
      
    }
    
    if(kernelY$method=="DTW"){
      
      ### default values for missing options ###
      
      if(is.null(kernelY$fam)) kernelY$fam <- "rbf"
      
      ### removal of undesired options in kernelY ###
      
      expected.kY <- c("method", "fam")
      found.kY <- names(kernelY)
      if(!(all(found.kY %in% expected.kY))){
        undesired.kY <- setdiff(found.kY, expected.kY)
        kernelY[undesired.kY] <- NULL
        warning("The list kernelY can only contain the following option: method. \n",
                "Other options have been removed.", call.=FALSE)
      }

      ### check of remaining options in kernelY ###
      
      # kernelY$fam must be a string.
      # the string in kernelY$fam must belong to L$dtw.kernels
      
      if(!is.character(kernelY$fam) | length(kernelY$fam)>1){
        stop("kernelY$fam must be a string.", call.=FALSE)
      }else{
        if(!(kernelY$fam %in% L$dtw.kernels)){
          stop("Invalid kernel name in kernelY$fam. \n", 
               "Available kernels include: ",
               paste0(L$dtw.kernels, collapse=', '), ".", call.=FALSE)
        }
      }

    }

  }
  
  ### For paramY ###############################################################
  
  ### minimum requirements ###
  
  # 1) paramY must be a scalar value or a vector of values
  # 2) elements in paramY must be NA or numeric
  
  msg <- "paramY must be a scalar value or a vector of values."
  
  if(is.numeric(paramY)){
    
    if(!is.null(dim(paramY))) stop(msg, call.=FALSE)
    
  }else if(is.logical(paramY)){
    
    if(!all(is.na(paramY))) stop(msg, call.=FALSE)
    
  }else{
    
    stop(msg, call.=FALSE)
    
  }

  ### For consistency between kernelY and paramY ###############################
  
  # --> see tell.sensiHSIC
  
  ### For estimator.type #######################################################

  # estimator.type must be either "V-stat" or "U-stat"
  
  if(!is.character(estimator.type) | length(estimator.type)>1){
    stop("estimator.type must be a string.", call.=FALSE)
  }else{
    if(!(estimator.type %in% c("V-stat", "U-stat"))){
      estimator.type <- "V-stat"
      warning("estimator.type must be: V-stat or U-stat. \n",
              "By default, estimator.type=V-stat has been selected.", call.=FALSE)
    }
  }
  
  ### For estimator.type and cond ##############################################
  
  # conditional R2-HSIC cannot be computed with U-statistics
  # only the conditional HSIC are thus estimated
  
  if(!is.null(cond) && estimator.type=="U-stat"){
    estimator.type <- "V-stat"
    warning("The conditional R2-HSIC indices cannot be estimated with U-statistics. \n",
            "By default, estimator.type=V-stat has been selected.", call.=FALSE)
  }
  
  ### For nboot ################################################################

  # nboot must be non-negative integer
  
  if(!is.numeric(nboot) | length(nboot)>1){
    stop("nboot must be a non-negative integer.", call.=FALSE)
  }else{
    if(nboot<0){
      stop("nboot must be a non-negative integer.", call.=FALSE)
    }else{
      if(!(nboot%%1==0)){
        nboot <- floor(nboot)
        warning("nboot was not an integer. \n",
                "floor(nboot) has been taken instead.", call.=FALSE)
      }
    }
  }
  
  ### For conf #################################################################
  
  # conf must be a real number in [0,1]
  
  if(!is.numeric(conf) | length(conf)>1){
    stop("conf must be a real number between 0 and 1.", call.=FALSE)
  }else{
    if(conf<0 | conf>1){
      conf <- 0.95
      warning("conf was not between 0 and 1. \n",
              "By default, conf=0.95 has been selected.", call.=FALSE)
    }
  }
  
  ### For anova ################################################################
    
  ### minimum requirements ###
  
  # 1) anova must be a list of options
  # 2) anova must contain at least one element named "obj"
  
  if(!is.list(anova)){
    stop("anova must be a list of options.", call.=FALSE)
  }else{
    if(is.null(anova$obj)){
      stop("Objective not defined. You must specify anova$obj.", call.=FALSE)
    }
  }
  
  ### default value for the other option ###
  
  if(is.null(anova$is.uniform)) anova$is.uniform <- TRUE
  
  ### removal of undesired options in anova ###
  
  expected.anova <- c("obj", "is.uniform")
  found.anova <- names(anova)
  if(!(all(found.anova %in% expected.anova))){
    undesired.anova <- setdiff(found.anova, expected.anova)
    anova[undesired.anova] <- NULL
    warning("The list anova can only contain the following options: ",
            paste0(expected.anova, collapse=', '), ". \n", 
            "Other options have been removed.", call.=FALSE)
  }
  
  ### check of remaining options in anova ###
  
  # anova$obj must be a string among "no", "FO", "TO" and "both"
  
  list.obj <- c("no", "FO", "TO", "both")
  if(!is.character(anova$obj) | length(anova$obj)>1){
    stop("anova$obj must be a string.", call.=FALSE)
  }else{
    if(!(anova$obj %in% list.obj)){
      anova$obj = "no"
      warning("anova$obj does not belong to: ",
              paste0(list.obj, collapse=', '), ". \n",
              "By default, anova$obj=no has been selected.", call.=FALSE)
    }
  }
  
  # anova$is.uniform must be boolean
  
  if(length(anova$is.uniform)>1){
    stop("anova$is.uniform must be a boolean.", call.=FALSE)
  }else{
    if(!(isTRUE(anova$is.uniform) | isFALSE(anova$is.uniform))){
      stop("anova$is.uniform must be a boolean.", call.=FALSE)
    }
  }
  
  # if ANOVA is required, all input kernels must be in L$anova.kernels or in L$conv.kernels
  
  if(!(anova$obj=="no")){
    for(i in 1:nkX){
      ker <- kernelX[i]
      if(!(ker %in% L$anova.kernels)){
        if(ker %in% L$conv.kernels){
          # first case: ker can be converted into an ANOVA kernel
          kernelX[i] <- paste0(ker, "_anova")
          paramX[i] <- NA
          warning(ker, " is not an ANOVA kernel. Its variant ", 
                  ker, "_anova has been taken instead.", call.=FALSE)
        }else{
          # second case: ker cannot be converted into an ANOVA kernel 
          kernelX[i] <- "sobolev1"
          paramX[i] <- NA
          warning(ker, " is not an ANOVA kernel and it cannot be converted into an ANOVA kernel. \n",
                  "By default, sobolev1 has been selected.", call.=FALSE)
          
        }
      }
    }
  }

  # the HSIC-ANOVA decomposition cannot be achieved for conditional SA
  # after conditioning, the input variables are non longer independent

  if(!is.null(cond) && anova$obj!="no"){
    anova <- list(obj="no", check=FALSE)
    warning("cond and anova cannot be enabled simultaneously. \n",
            "anova has been ignored.", call.=FALSE)
  }
  
  ### For sensi ################################################################
  
  if(!is.null(sensi)){
    
    # sensi must be an object of class sensiHSIC
    # sensi should contain an object named KX or an object named KY 
    
#    if(!(class(sensi)=="sensiHSIC")){ # bug
    if(!(inherits(sensi,"sensiHSIC"))){
      stop("sensi must be an object of class sensiHSIC.", call.=FALSE)
    }else{
      if(is.null(sensi$KX) && is.null(sensi$KY)){
        warning("The sensi option is useless here. ",
                "Neither KX nor KY is stored in the provided object. \n",
                "All Gram matrices are reassembled from scratch.", call.=FALSE)
      }
    }
    
  }

  ### For save.GM ##############################################################
  
  ### requirements ###
  
  # 1) save.GM must be a list of options
  # 2) save.GM must contain one element named "KX" and one element named "KY"
  # 3) save.GM$KX and save.GM$KY must be boolean
  
  # must be list of options
  if(!is.list(save.GM)){
    stop("save.GM must be a list of options.", call.=FALSE)
  }else{
    # must contain one element named KX
    if(is.null(save.GM$KX)){
      stop("Storage of KX not decided. save.GM$KX must be specified.", call.=FALSE)
    }else{
      # must be of length 1
      if(length(save.GM$KX)>1){
        stop("save.GM$KX must be a boolean.", call.=FALSE)
      }else{
        # must be either TRUE or FALSE
        if(!(isTRUE(save.GM$KX) | isFALSE(save.GM$KX))){
          stop("save.GM$KX must be a boolean.", call.=FALSE)
        }
      }
    }
    # must contain one element named KY
    if(is.null(save.GM$KY)){
      stop("Storage of KY not decided. save.GM$KY must be specified.", call.=FALSE)
    }else{
      # must be of length 1
      if(length(save.GM$KY)>1){
        stop("save.GM$KY must be a boolean.", call.=FALSE)
      }else{
        # must be either TRUE or FALSE
        if(!(isTRUE(save.GM$KY) | isFALSE(save.GM$KY))){
          stop("save.GM$KY must be a boolean.", call.=FALSE)
        }
      }
    }
    ### removal of undesired options in save.GM ###
    expected.GM <- c("KX", "KY")
    found.GM <- names(save.GM)
    if(!(all(found.GM %in% expected.GM))){
      undesired.GM <- setdiff(found.GM, expected.GM)
      save.GM[undesired.GM] <- NULL
      warning("The list save.GM can only contain the following options: ",
              paste0(expected.GM, collapse=', '), ". \n", 
              "Other options have been removed.", call.=FALSE)
    }
  }
  
  #########################
  ### Main instructions ###
  #########################
  
  x <- list(model=model, X=X, target=target, cond=cond, 
           kernelX=kernelX, paramX=paramX,
           kernelY=kernelY, paramY=paramY,
           estimator.type=estimator.type,
           nboot=nboot, conf=conf,
           anova=anova,
           sensi=sensi, save.GM=save.GM,
           call = match.call())

  class(x) <- "sensiHSIC"
  
  if(!is.null(x$model)){ 
    response(x, other_types_allowed=TRUE, ...)
    tell(x, ...)
  }

  return(x)
  
}

### tell.sensiHSIC #############################################################

tell.sensiHSIC <- function(x, y = NULL, ...){
  
  id <- deparse(substitute(x))

  # How to start from given data?
  
  # dimensions
  n <- nrow(x$X)
  p <- ncol(x$X)
  
  ### useful lists #############################################################
  
  L <- available.material()
  
  ### For y ####################################################################
  
  # 1) if y is not provided: 
  # --> a) it must be checked that x$y is available
  # --> b) the output matrix Y = x$y is defined
  # 2) if y is provided:
  # --> a) y must be a matrix or a data frame
  # --> b) y must contain as many observations as X
  # --> c) the output matrix Y = y is defined
  
  if(is.null(y)){ # y is not provided
    if(is.null(x$y)) stop("y not found.", call.=FALSE)
    # if y is a vector --> y is converted into a n-by-1 matrix
    if(is.numeric(x$y) && is.null(dim(x$y))){
      Y <- as.matrix(x$y)
    }else{
      Y <- x$y
    }
  }else{ # y is provided
    # y must be a matrix or must be converted into a matrix:
    # if y is a vector --> y is converted into a n-by-1 matrix
    # if y is a data frame --> y is converted into a matrix
    if(is.numeric(y) && is.null(dim(y))){ # y is a vector
      y <- as.matrix(y)
    }else if(is.data.frame(y)){ # y is a data.frame
      y <- as.matrix(unname(y))
    }else if(!is.matrix(y)){
      stop("y must be a matrix or a data frame.", call.=FALSE)
    }
    # must contain as many observations as X
    if(!(nrow(y)==n)) stop("y must contain as many observations as X.", call.=FALSE)
    Y <- y
  }
  
  q <- ncol(Y)
  
  ### For consistency between q, target, cond and kernelY ######################
  
  # Target SA is only possible if the output is scalar
  if(q>1 && !is.null(x$target)){
    stop("Target SA is only possible if the output is scalar.", call.=FALSE)
  }
  
  # Conditional SA is only possible if the output is scalar
  if(q>1 && !is.null(x$cond)){
    stop("Conditional SA is only possible if the output is scalar.", call.=FALSE)
  }
  
  # Methods based on PCA, GAK or PCA are specific to non-scalar outputs
  if(q==1 && is.list(x$kernelY)){
    
    msg <- switch(x$kernelY$method,
                  "PCA"={ "Preliminary dimension reduction based on PCA" },
                  "GAK"={ "Using the global aligmnent kernel (GAK)" },
                  "DTW"={ "Using dynamic time warping (DTW)" },
                  NA)
    
    msg <- paste0(msg, " is not appropriate for a scalar output. \n",
                 "kernelY must be a string specifying the selected kernel. \n",
                 "Available kernels include: ", 
                 paste0(L$kernels, collapse=', '), ".")
    
    stop(msg, call.=FALSE)
    
  } 

  # How to redefine x$y?
  # --> for basic SA or for CSA, x$y is simply Y
  # --> for TSA, x$y must be computed with the function weightTSA
  if(is.null(x$target)){
    x$y <- Y
  }else{
    x$y <- matrix(weightTSA(Y, x$target$c, x$target$upper, x$target$type, x$target$param))
  }
  
  ### For target ###############################################################
  
  ### additional check ###
  
  # let us imagine that target is enabled
  # --> if the weight function is "indicTh": kernelY must be "categ"
  # --> if the weight function is "zeroTh", "logistic" or "exp1side": kernelY must not be "categ"
  
  msg <- paste0("Be careful. The output kernel and the weight function are not consistent. \n",
                "If target$type=indicTh, kernelY must refer to the categorical kernel. \n",
                "In all other cases, kernelY must refer to a continuous kernel. \n")
  
  if(!is.null(x$target)){
    if(x$target$type=="indicTh"){
      if(!(x$kernelY=="categ")){
        x$kernelY <- "categ"
        warning(msg, "By default, kernelY=categ has been selected.", call.=FALSE)
      }
    }else{
      if(x$kernelY=="categ"){
        x$kernelY <- "rbf"
        warning(msg, "By default, kernelY=rbf has been selected.", call.=FALSE)
      }
    }
  }

  ### For anova ################################################################

  # only if the HSIC-ANOVA decomposition is required
  
  if(x$anova$obj!="no"){

    if(x$anova$is.uniform==TRUE){
    
      # if anova$is.uniform=TRUE, X is examined 
      # if some samples do not lie in [0,1] --> error

      if(!(all(x$X>=0 & x$X<=1))){
        stop("You wrongly specified is.uniform=TRUE in anova. \n",
             "In fact, some samples in X do not lie in [0,1]. You can either: \n",
             "1) Rescale the data in X so that all samples lie in [0,1]. \n",
             "Be careful: this may lead to change the function given to model. \n",
             "2) Specify is.uniform=FALSE in anova.", call.=FALSE)
      }
      
    }else{
      
      # if anova$is.uniform=FALSE, non-parametric rescaling is operated
      
      x$X <- apply(x$X, 2, function(col){ rank(col) })/n
      
    }
    
  }

  ### For kernelX and paramX ###################################################
  
  ### conversion into fully specified objects ###
  
  convX <- kernel.param.X(x$kernelX, x$paramX, x$X)
  x$kernelX <- convX$kernelX
  x$paramX <- convX$paramX
  
  ### For kernelY and paramY ###################################################
  
  ### additional checks based on q ###
  
  # 1) if kernelY is a string or a vector of strings:
  # --> kernelY must be of length either 1 or q
  # --> paramY must be of length either 1 or q
  
  # 2) if kernelY is a list of options:
  # a) if kernelY$method="PCA"
  # --> kernelY$fam must be of length 1
  # --> paramY must be NA
  # b) if kernelY$method="GAK"
  # --> paramY must be of length 2
  # c) if kernelY$method="DTW"
  # --> kernelY$fam must be of length 1
  # --> paramY must be NA
  
  npY <- length(x$paramY)
  
  if(is.character(x$kernelY)){
    
  # kernelY is a string or a vector of strings
    
    nkY <- length(x$kernelY)
    # check of kernelY
    if(!(nkY==1 | nkY==q)){
      stop("kernelY must be of length 1 or q (number of output variables).", call.=FALSE)
    }
    # check of paramY
    if(!(npY==1 | npY==q)){
      stop("paramY must be of length 1 or q (number of output variables).", call.=FALSE)
    }
    
  }else{
    
  # kernelY is a list of options
    
    # case kernelY$method="PCA"
    if(x$kernelY$method=="PCA"){
      # check of kernelY
      nkY <- length(x$kernelY$fam)
      if(!(nkY==1)){
        stop("If kernelY$method=PCA, kernelY$fam must be of length 1.", call.=FALSE)
      }
      # check of paramY
      if(!(npY==1)){
        stop("If kernelY$method=PCA, paramY must be NA.", call.=FALSE)
      }else{
        if(!(is.na(x$paramY))){
          x$paramY <- NA
          warning("If kernelY$method=PCA, paramY is computed automatically. \n",
                  "Default assignment to NA has been accepted.", call.=FALSE)
        }
      }
    }

    # case kernelY$method="GAK"
    if(x$kernelY$method=="GAK"){
      # check of paramY
      if(npY==1){
        if(is.na(x$paramY)){
          x$paramY <- rep(NA, 2)
        }else{
          stop("If kernelY$method=GAK, paramY must be of length 2.", call.=FALSE)
        }
      }else if(npY>=3){
        stop("If kernelY$method=GAK, paramY must be of length 2.", call.=FALSE)
      }
    }
    
    # case kernelY$method="DTW"
    if(x$kernelY$method=="DTW"){
      # check of kernelY
      nkY <- length(x$kernelY$fam)
      if(!(nkY==1)){
        stop("If kernelY$method=DTW, kernelY$fam must be of length 1.", call.=FALSE)
      }
      # check of paramY
      if(!(npY==1)){
        stop("If kernelY$method=DTW, paramY must be NA.", call.=FALSE)
      }else{
        if(!(is.na(x$paramY))){
          x$paramY <- NA
          warning("If kernelY$method=DTW, paramY must be NA.", call.=FALSE)
        }
      }
    }
    
  }
  
  ### conversion into fully specified objects ###
  
  convY <- kernel.param.Y(x$kernelY, x$paramY, x$y)
  x$kernelY <- convY$kernelY
  x$paramY <- convY$paramY
  
  #########################
  ### Main instructions ###
  #########################
  
  ### computation of Gram matrices ###
  
  # if an object sensi is provided:
  # a) if KX is stored:
  # --> the input Gram matrices must not be recomputed
  # --> kernelX, paramX and KX must be taken from sensi
  # b) if KY is stored:
  # --> the output Gram matrix must not be recomputed
  # --> kernelY, paramY and KY must be taken from sensi
  
  # --> sensi is eliminated at the end of the process
  
  need2cpt.KX <- TRUE
  need2cpt.KY <- TRUE
  
  if(!is.null(x$sensi)){
    # access at sensi$KX
    if(!is.null(x$sensi$KX)){
      need2cpt.KX <- FALSE
    }
    # access at sensi$KY
    if(!is.null(x$sensi$KY)){
      need2cpt.KY <- FALSE
    }
  }
  
  # computation of KX (if required)
  if(need2cpt.KX){
    KX <- compute.KX(x$X, x$kernelX, x$paramX)
  }else{
    x$kernelX <- x$sensi$kernelX
    x$paramX <- x$sensi$paramX
    KX <- x$sensi$KX 
    if(p==1){
      warning("The input Gram matrix has been extracted ", 
              "from an older object of class sensiHSIC.", call.=FALSE)
    }else{
      warning("The input Gram matrices have been extracted ", 
              "from an older object of class sensiHSIC.", call.=FALSE)
    }
  }
  
  # computation of KY (if required)
  if(need2cpt.KY){
    GM.Y <- compute.KY(x$y, x$kernelY, x$paramY)
    x$kernelY <- GM.Y$kernelY
    x$paramY <- GM.Y$paramY
    KY <- GM.Y$KY
  }else{
    x$kernelY <- x$sensi$kernelY
    x$paramY <- x$sensi$paramY
    KY <- x$sensi$KY
    warning("The output Gram matrix has been extracted ", 
            "from an older object of class sensiHSIC.", call.=FALSE)
  }
  
  x$sensi <- NULL

  ### storage of Gram matrices ###
  
  if(x$save.GM$KX){
    x$KX <- KX
  }
  
  if(x$save.GM$KY){
    x$KY <- KY
  }

  ### required statistics according to the objective ###
  
  stats <- switch(x$anova$obj, 
                  "no"={ c("HSICXY", "S") },
                  "FO"={ c("HSICXY", "S", "FO", "denom") }, 
                  "TO"={ c("HSICXY", "S", "TO", "TO.num", "denom") },
                  "both"={ c("HSICXY", "S", "FO", "TO", "TO.num", "denom") },
                  NA)
  
  ### computation of HSIC indices ###

  data <- cbind(x$X, x$y)
  
  if(!is.null(x$cond)){
    
    w <- weightTSA(x$y, x$cond$c, x$cond$upper, x$cond$type, x$cond$param)
    # if w the zero vector, the conditioning operation is useless
    if(all(w==0)){
      stop("The conditioning event is never observed among the provided dataset. \n",
           "Please check if cond is well-defined.", call.=FALSE)
    }
    x$weights <- w/mean(w)
    
  }
  
  if(x$nboot==0){

    if(is.null(x$cond)){
      res <- estim.sensiHSIC(KX, KY, x$estimator.type, stats)
    }else{
      res <- estim.cond.sensiHSIC(KX, KY, x$weights, stats)
    }

    x$HSICXY <- data.frame(res$HSICXY)
    rownames(x$HSICXY) <- paste("X", 1:p, sep="")
    colnames(x$HSICXY) <- "original"
    
    if(is.null(x$cond) | x$estimator.type=="V-stat"){
      x$S <- data.frame(res$S)
      rownames(x$S) <- paste("X", 1:p, sep="")
      colnames(x$S) <- "original"
    }
    
    if(x$anova$obj %in% c("FO", "TO", "both")){

      if(x$anova$obj!="TO"){
        x$FO <-  data.frame(res$FO)
        rownames(x$FO) <- paste("X", 1:p, sep="")
        colnames(x$FO) <- "original"
      }
      
      if(x$anova$obj!="FO"){
        x$TO <-  data.frame(res$TO)
        x$TO.num <-  data.frame(res$TO.num)
        rownames(x$TO) <- rownames(x$TO.num) <- paste("X", 1:p, sep="")
        colnames(x$TO) <- colnames(x$TO.num) <- "original"
      }
      
      x$denom <-  data.frame(res$denom)
      colnames(x$denom) <- "original"
      
    }

  }else{
    
    if(is.null(x$cond)){
      
      all.boot <- boot(data, statistic=function(data, ind){ 
        estim.sensiHSIC(KX[ind, ind, ], KY[ind, ind], 
                        x$estimator.type, stats, vector=TRUE) 
        }, R=x$nboot)
      
    }else{
      
      all.boot <- boot(data, statistic=function(data, ind){ 
        estim.cond.sensiHSIC(KX[ind, ind, ], KY[ind, ind], x$weights[ind], stats, vector=TRUE) 
        }, R=x$nboot)
      
    }
    
    res <- bootstats(all.boot, x$conf, "basic")
    
    x$HSICXY <- res[1:p, ] # HSIC indices
    rownames(x$HSICXY) <- paste("X", 1:p, sep="")

    if(is.null(x$cond) | x$estimator.type=="V-stat"){
      x$S <- res[(p+1):(2*p),] # R2-HSIC indices
      rownames(x$S) <- paste("X", 1:p, sep="")
    }

    if(x$anova$obj=="FO"){
      x$FO <- res[(2*p+1):(3*p), ]      # first-order HSIC indices
      x$denom <- res[3*p+1, ]           # common denominator
      rownames(x$FO) <- paste("X", 1:p, sep="")
      rownames(x$denom)  <- ""
    }
    
    if(x$anova$obj=="TO"){
      x$TO <- res[(2*p+1):(3*p), ]      # total-order HSIC indices
      x$TO.num <- res[(3*p+1):(4*p), ]  # numerators of total-order HSIC indices
      x$denom <- res[4*p+1, ]           # common denominator
      rownames(x$TO) <- rownames(x$TO.num) <- paste("X", 1:p, sep="")
      rownames(x$denom)  <- ""
    }
    
    if(x$anova$obj=="both"){
      x$FO <- res[(2*p+1):(3*p), ]      # first-order HSIC indices
      x$TO <- res[(3*p+1):(4*p), ]      # total-order HSIC indices
      x$TO.num <- res[(4*p+1):(5*p), ]  # total-order HSIC indices
      x$denom <- res[5*p+1, ]           # common denominator
      rownames(x$FO) <- rownames(x$TO) <- rownames(x$TO.num) <- paste("X", 1:p, sep="")
      rownames(x$denom)  <- ""
    }
    
  }
  
  assign(id, x, parent.frame())
}

### estim.sensiHSIC ############################################################

estim.sensiHSIC <- function(KX, KY, estimator.type, stats, vector=FALSE){

  ### inputs ###

  # KX:               n x n x p   array     (concatenation of input Gram matrices)
  # KY:               n x n       matrix    (output Gram matrix)
  # estimator.type:   1 x 1       character ("U-stat" or "V-stat")
  # stats:            1 x s       character (all desired statistics)
  # vector:           1 x 1       logical   (FALSE if the output must be a list)
  #                                         (TRUE ________________________ vector)

  # In stats, the strings must be ordered as follows:
  # "HSICXY", "S", "FO", "TO", "TO.num", "denom".

  ### output ###

  # 1) if vector=FALSE, res is a list of results (depending on stats)

  # HSICXY:     1 x p   numeric   (HSIC indices)
  # S:          1 x p   numeric   (R2-HSIC indices)
  # FO:         1 x p   numeric   (first-order HSIC-ANOVA indices)
  # TO:         1 x p   numeric   (total-order HSIC-ANOVA indices)
  # TO.num:     1 x p   numeric   (numerators of total-order HSIC-ANOVA indices)
  # denom:      1 x 1   numeric   (common denominator of HSIC-ANOVA indices)

  # 2) if vector=TRUE, res is a numeric of length L depending on stats

  ### For stats ################################################################

  if(vector){

    reference <- c("HSICXY", "S", "FO", "TO", "TO.num", "denom")

    s <- length(stats)
    order <- match(reference, stats)
    order <- order[!is.na(order)]

    if(!identical(order, 1:s)){
      stats <- stats[order]
      warning("The required statistics were sorted in the following order: \n",
              paste0(stats, collapse=', '), ".", call.=FALSE)
    }

  }

  ### Body of the function #####################################################

  n <- dim(KX)[1]
  p <- dim(KX)[3]

  res <- list() # pre-allocation of the output object

  ### computations that occur in several cases ###

  if(any(c("HSICXY", "S", "FO") %in% stats)){

    # estimation of all HSIC(Xi,Y)
    HSIC.XY <- rep(NA, p)
    for(i in 1:p){
      HSIC.XY[i] <- estim.HSIC(KX, KY, i, estimator.type)
    }

  }

  if(any(c("FO", "TO", "TO.num", "denom") %in% stats)){

    # estimation of the common denominator
    HSIC.denom <- estim.HSIC(KX, KY, 1:p, estimator.type)

  }

  if(any(c("TO", "TO.num") %in% stats)){

    # estimation of all HSIC(X_{1:d \ i},Y)
    HSIC.compl <- rep(NA, p)
    for(i in 1:p){
      HSIC.compl[i] <- estim.HSIC(KX, KY, setdiff(1:p, i), estimator.type)
    }

  }

  ### all other computations (and storage) ###

  if("HSICXY" %in% stats){
    res$HSICXY <- HSIC.XY # storage
  }

  if("S" %in% stats){

    # estimation of all HSIC(Xi,Xi)
    HSIC.XX <- rep(NA, p)
    for(i in 1:p){
      HSIC.XX[i] <- estim.HSIC(KX, KX[,,i], i, estimator.type)
    }

    # estimation of HSIC(Y,Y)
    LY <- replicate(1, KY, simplify="array")
    HSIC.YY <- estim.HSIC(LY, KY, 1, estimator.type)

    # estimation of R2 HSIC indices
    R2.HSIC <- HSIC.XY/sqrt(HSIC.YY*HSIC.XX)

    res$S <- R2.HSIC # storage

  }

  if("FO" %in% stats){

    # estimation of first-order HSIC-ANOVA indices
    HSIC.FO <- HSIC.XY/HSIC.denom

    res$FO <- HSIC.FO # storage

  }

  if("TO" %in% stats){

    # estimation of total-order HSIC-ANOVA indices
    HSIC.TO <- 1-HSIC.compl/HSIC.denom

    res$TO <- HSIC.TO # storage

  }

  if("TO.num" %in% stats){

    # estimation of the numerators of total-order HSIC-ANOVA indices
    HSIC.TO.num <- HSIC.denom-HSIC.compl

    res$TO.num <- HSIC.TO.num # storage

  }

  if("denom" %in% stats){

    res$denom <- HSIC.denom # storage

  }

  ### conversion from a list into a vector (if required) ###

  if(vector){
    res <- unlist(res, use.names=FALSE)
  }

  return(res)

}

### estim.HSIC #################################################################
    
estim.HSIC <- function(KX, KY, A, estimator.type){
  
  ### inputs ###
  
  # KX:               n x n x p   array     (concatenation of input Gram matrices)
  # KY:               n x n       matrix    (output Gram matrix)
  # A:                1 x r       numeric   (subset of 1,...,p)
  # estimator.type:   1 x 1       character ("U-stat" or "V-stat")
  
  ### output ###
  
  # HSIC              1 x 1       numeric   (HSIC index estimator)
  
  n <- dim(KX)[1]
  p <- dim(KX)[1]
  
  PX <- 1
  
  for(i in A){
    PX <- PX*KX[,,i]
  }
  
  HSIC <- HSIC.fun(PX, KY, estimator.type)
  
  return(HSIC)
  
}

### HSIC.fun ###################################################################

HSIC.fun <- function(KX, KY, estimator.type){
  
  ### inputs ###
  
  # KX:               n x n       matrix    (input Gram matrix)
  # KY:               n x n       matrix    (output Gram matrix)
  # estimator.type:   1 x 1       character ("U-stat" or "V-stat")
  
  ### output ###
  
  # HSIC              1 x 1       numeric   (HSIC index estimator)
  
  n <- dim(KX)[1]
  
  if(estimator.type=="V-stat"){ # formula for the V-statistic
    
    UXY <- sum(KX*KY)
    VXY <- sum(colSums(KX)%*%KY)
    WXY <- sum(KX)*sum(KY)
    HSIC <- UXY/(n^2) - 2*VXY/(n^3) + WXY/(n^4)
    
  }else{ # formula for the U-statistic
    
    diag(KX) <- diag(KY) <- 0
    
    UXY <- sum(KX*KY)
    VXY <- sum(colSums(KX)%*%KY) - UXY
    WXY <- sum(KX)*sum(KY) - 2*UXY - 4*VXY
    HSIC <- UXY/(n*(n-1)) - 2*VXY/(n*(n-1)*(n-2)) + WXY/(n*(n-1)*(n-2)*(n-3))
    
  }
  
  return(HSIC)
  
}

### estim.cond.sensiHSIC #######################################################

estim.cond.sensiHSIC <- function(KX, KY, weights, stats, vector=FALSE){

  ### inputs ###
  
  # KX:               n x n x p   array     (concatenation of input Gram matrices)
  # KY:               n x n       matrix    (output Gram matrix)
  # weights:          1 x n       numeric   (conditioning weights after normalization)
  # stats:            1 x s       character (all desired statistics)
  # vector:           1 x 1       logical   (FALSE if the output must be a list)
  #                                         (TRUE ________________________ vector)
  
  # In stats, the strings must be ordered as follows: "HSICXY", "S".
  
  ### output ###
  
  # 1) if vector=FALSE, res is a list of results (depending on stats)

  # HSICXY:     1 x p   numeric   (HSIC indices)
  # S:          1 x p   numeric   (R2-HSIC indices)

  # 2) if vector=TRUE, res is a numeric of length L depending on stats
  
  ### For stats ################################################################
  
  if(vector){
    
    reference <- c("HSICXY", "S")
    if(identical(stats, rev(reference))){
      stats <- rev(stats)
      warning("The required statistics were sorted in the following order: ", 
              "HSICXY, S.", call.=FALSE)
    }

  }
  
  ### Body of the function #####################################################

  n <- dim(KX)[1]
  p <- dim(KX)[3]
  
  res <- list()
  
  ### computations that occur in several cases ###
  
  # estimation of all HSIC(Xi,Y)
  HSIC.XY <- rep(NA, p)
  for(i in 1:p){
    HSIC.XY[i] <- cond.HSIC.fun(KX[,,i], KY, weights)
  }

  ### all other computations (and storage) ###
  
  if("HSICXY" %in% stats){
    res$HSICXY <- HSIC.XY # storage
  }
  
  if("S" %in% stats){
    
    # estimation of all HSIC(Xi,Xi)
    HSIC.XX <- rep(NA, p)
    for(i in 1:p){
      HSIC.XX[i] <- cond.HSIC.fun(KX[,,i], KX[,,i], weights)
    }
    
    # estimation of HSIC(Y,Y)
    HSIC.YY <- cond.HSIC.fun(KY, KY, weights)
    
    # estimation of R2 HSIC indices
    R2.HSIC <- HSIC.XY/sqrt(HSIC.YY*HSIC.XX)
    
    res$S <- R2.HSIC # storage
    
  }
  
  ### conversion from a list into a vector (if required) ###
  
  if(vector){
    res <- unlist(res, use.names=FALSE)
  }
  
  return(res)

}

### cond.HSIC.fun ##############################################################

cond.HSIC.fun <- function(KX, KY, weights){
  
  ### inputs ###
  
  # KX:               n x n   matrix      (concatenation of input Gram matrices)
  # KY:               n x n   matrix      (output Gram matrix)
  # weights:          1 x n   numeric     (conditioning weights after normalization)
  
  ### output ###
  
  # C.HSIC              1 x 1       numeric   (C-HSIC index estimator)
  
  n <- dim(KX)[1]
  
  wu <- matrix(rep(weights, n), n, n) 
  W <- wu*t(wu)
  # Rk: wu corresponds to the matrix product WU
  # --> see the paper by Marrel & Chabridon (2021)
  
  # formula for the V-statistic
  UXY <- sum(KX*KY*W)
  VXY <- sum((colSums(KX*wu)%*%(KY*W)))
  WXY <- sum(KX*W)*sum(KY*W)
  
  C.HSIC <- UXY/(n^2) - 2*VXY/(n^3) + WXY/(n^4)

  return(C.HSIC)
  
}

### check.param ################################################################

check.param <- function(ker, par, X){
  
  ### inputs ###
  
  # ker:      1 x 1   character     (kernel family)
  # par:      1 x 1   numeric/NA    (kernel parameter)
  # X:        1 x n   numeric       (samples)
  
  ### output ###
  
  # par:      1 x 1   numeric     (kernel parameter after checking)
  
  ### useful lists ###
  
  L <- available.material()
  
  ### rbf-like kernels ###
  
  if(ker %in% L$sd.kernels){
    if(is.na(par)){
      par <- sd(X)
    }else{
      if(par<=0){
        par <- sd(X)
        warning("For the ", ker, " kernel, the parameter must be a positive real number. \n",
                "Default assignment has been accepted.", call.=FALSE)
      }
    }
  }
  
  ### kernels with no parameter ###
  
  if(ker %in% L$free.kernels){
    if(!is.na(par)){
      par <- NA
      warning("The ", ker, " kernel does not have parameter. \n",
              "Default assignment has been accepted.", call.=FALSE)
    }
  }
  
  ### dcov kernel ###
  
  if(ker=="dcov"){
    if(is.na(par)){
      par <- 1
    }else{
      if(par<=0){
        par <- 1
        warning("For the dcov kernel, the parameter must be a positive real number. \n",
                "Default assignment has been accepted.", call.=FALSE)
      }
    }
  }
  
  ### categorical kernel ###
  
  if(ker=="categ"){
    if(is.na(par)){
      par <- 0
    }else{
      if(!(par %in% c(0,1))){
        par <- 0
        warning("For the categ kernel, the parameter must be 0 or 1. \n",
                "Default assignment has been accepted.", call.=FALSE)
      }
    }
  }
  
  return(par)
  
}

### kernel.param.X #############################################################

kernel.param.X <- function(kernelX, paramX, X){
  
  ### inputs ###
  
  # kernelX:    1 x 1   character     (one kernel family)
  #             1 x p   character     (kernel families)
  
  # paramX:     1 x 1   numeric/NA    (one parameter value)
  #             1 x p   numeric/NA    (parameter values)
  
  # X:          n x p   matrix        (input samples)
  
  ### outputs ###
  
  # kernelX:    1 x p   character     (kernel families)
  # paramX:     1 x p   numeric/NA    (parameter values)
  
  p <- ncol(X)            # nb of input variables
  nkX <- length(kernelX)  # nb of specified kernels
  npX <- length(paramX)   # nb of specified parameters
  
  if(nkX==1) kernelX <- rep(kernelX, p)
  if(npX==1) paramX <- rep(paramX, p)
  
  for(i in 1:p){
    paramX[i] <- check.param(kernelX[i], paramX[i], X[,i])
  }
  
  conv <- list(kernelX=kernelX, paramX=paramX)
  return(conv)
  
}

### kernel.param.Y #############################################################

kernel.param.Y <- function(kernelY, paramY, Y){
  
  ### inputs ###
  
  # kernelY:    1 x 1     character     (one kernel family)
  #             1 x q     character     (kernel families)
  #             1 x ...   list          (options describing the kernel)
  
  # paramY:     1 x 1     numeric/NA    (one parameter value)
  #             1 x q     numeric/NA    (parameter values)
  #             1 x 2     numeric/NA    (parameter values for the GA kernel)
  
  # Y:          n x q     matrix        (output samples)
  
  ### outputs ###
  
  # kernelY:    1 x q     character     (kernel families)
  #             1 x ...   list          (options describing the kernel)   

  # paramY:     1 x q     numeric/NA    (parameter values)
  #             1 x 2     numeric/NA    (parameter values for the GA kernel)         
  
  q <- ncol(Y)
  
  ### case where kernelY is a string or a vector of strings ###
  
  if(is.character(kernelY)){
    
    nkY <- length(kernelY)  # nb of specified kernels
    npY <- length(paramY)   # nb of specified parameters
    
    if(nkY==1) kernelY <- rep(kernelY, q)
    if(npY==1) paramY <- rep(paramY, q)

    for(i in 1:q){
      paramY[i] <- check.param(kernelY[i], paramY[i], Y[,i])
    }

  }
  
  ### case where kernelY is a list of options ###
  
  # if kernelY$method="PCA", nothing has to be done at this stage
  # if kernelY$method="GAK", check of the two parameters
  # if kernelY$method="DTW", there is nothing to do
  
  if(is.list(kernelY)){

    if(kernelY$method=="GAK"){
      
      # check of paramY[1]
      if(is.na(paramY[1])){
        paramY[1] <- median(dist(Y))*sqrt(q)
      }else{
        if(paramY[1]<=0){
          paramY[1] <- median(dist(Y))*sqrt(q)
          warning("For the global alignement kernel (GAK), ", 
                  "the first parameter must be a positive real number. \n",
                  "Default assignment has been accepted.", call.=FALSE)
        }
      }
      # check of paramY[2]
      if(is.na(paramY[2])){
        paramY[2] <- 0.2*q
      }else{
        if(paramY[2]<=0){
          paramY[2] <- 0.2*q
          warning("For the global alignement kernel (GAK), ", 
                  "the second parameter must be a positive real number. \n",
                  "Default assignment has been accepted.", call.=FALSE)
        }
      }
      
    }
    
  }
  
  conv <- list(kernelY=kernelY, paramY=paramY)
  return(conv)
  
}

### compute.KX #################################################################

compute.KX <- function(X, kernelX, paramX){
  
  ### inputs ###
  
  # X:          n x p     matrix        (input samples)
  # kernelX:    1 x p     character     (kernel families)
  # paramX:     1 x p     numeric/NA    (kernel parameters)
  
  ### output ###
  
  # KX:         n x n x d   array     (concatenation of input Gram matrices)
  
  n <- nrow(X)
  p <- ncol(X)
  
  KX <- array(NA, c(n,n,p))
  
  for(i in 1:p){
    KX[,,i] <- do.call(get(paste(kernelX[i], "_hsic", sep="")), 
                       list(x=X[,i], param=paramX[i]))
  }
  
  return(KX)
  
}

### compute.KY #################################################################

compute.KY <- function(Y, kernelY, paramY){

  ### inputs ###

  # Y:          n x q     matrix        (output samples)
  
  # kernelY:    1 x q     character     (kernel families)
  #             1 x ...   list          (options describing the kernel)
  
  # paramY:     1 x q     numeric/NA    (parameter values)
  #             1 x 2     numeric/NA    (parameter values for the GA kernel)

  ### output ###

  # KY:         n x n     matrix        (output Gram matrix)

  n <- nrow(Y)
  q <- ncol(Y)
  
  if(is.character(kernelY)){
    
    KY <- 1
    for(i in 1:q){
      KY <- KY*do.call(get(paste(kernelY[i], "_hsic", sep="")),
                       list(x=Y[,i], param=paramY[i]))
    }

  }else{
    
    if(kernelY$method=="PCA"){

      res.PCA <- PCA.KY(Y, kernelY)
      KY <- res.PCA$KY
      kernelY <- res.PCA$kernelY
      paramY <- res.PCA$paramY
      
    }else if(kernelY$method=="GAK"){
      
      KY <- GAK.KY(Y, paramY)
      
    }else if(kernelY$method=="DTW"){
      
      KY <- DTW.KY(Y, kernelY)
        
    }
    
  }

  res <- list(KY=KY, kernelY=kernelY, paramY=paramY)

  return(res)

}

### PCA.KY #####################################################################

PCA.KY <- function(Y, kernelY){
  
  ### inputs ###
  
  # Y:          n x q     matrix        (output samples)
  # kernelY:    1 x ...   list          (options describing the kernel)
  
  ### outputs ###
  
  # KY:         n x n     matrix        (output Gram matrix)
  # kernelY:    1 x ...   list          (options describing the kernel)  
  # paramY:     1 x ...   numeric       (parameters of the output kernel)
  
  n <- dim(Y)[1]
  
  # constant outputs must be removed
  mask.cst <- apply(Y, 2, sd, na.rm=TRUE)==0
  Y <- Y[,!mask.cst]
  
  # principal component analysis (PCA)
  pca.res <- prcomp(Y, retx=TRUE, 
                    center=kernelY$data.centering, # whether data need to be centered
                    scale.=kernelY$data.scaling)   # whether data need to be scaled
  lambda <- pca.res$sdev^2    # eigenvalues
  V <- unname(pca.res$x)      # eigenvectors
  
  if(!is.na(kernelY$expl.var)){ # expl.var is provided (but PC is NULL)
    
    # increasing sequence of variances
    cum.var <- cumsum(lambda)/sum(lambda)
    # nb of components to reach the desired level of output variance
    kernelY$PC <- min(which(cum.var>kernelY$expl.var))
    # percentage of explained variance
    kernelY$expl.var <- cum.var[kernelY$PC]
    
  }else{ # PC is provided (expl.var is NULL)
    
    # percentage of explained variance
    kernelY$expl.var <- sum(lambda[1:kernelY$PC])/sum(lambda)
    
  }
  
  PC <- kernelY$PC
  U <- V[, 1:PC, drop=FALSE]
  if(kernelY$position %in% c("intern", "extern")){
    ratios <- lambda[1:PC]/sum(lambda[1:PC])
  }
  
  # default assignment of kernel parameters
  param <- rep(NA, PC)
  for(i in 1:PC){
    param[i] <- check.param(kernelY$fam, NA, U[,i])
  }
  
  # computation of the final Gram matrix
  # 6 different cases:
  # a) kernelY$position="nowhere" + kernelY$combi="sum"
  # b) kernelY$position="nowhere" + kernelY$combi="prod"
  # c) kernelY$position="intern" + kernelY$combi="sum"
  # d) kernelY$position="intern" + kernelY$combi="prod"
  # e) kernelY$position="extern" + kernelY$combi="sum"
  # f) kernelY$position="extern" + kernelY$combi="prod" (not allowed)
  
  KU <- array(NA, c(n,n,PC))
  
  if(kernelY$position=="intern"){
    
    # computation of Gram matrices of the principal components
    # --> the weights are directly applied on the principal components
    
    for(i in 1:PC){
      KU[,,i] <- do.call(get(paste(kernelY$fam, "_hsic", sep="")),
                         list(x=ratios[i]*U[,i], param=param[i]))
    }
    
  }else if(kernelY$position=="extern"){
    
    # computation of Gram matrices of the principal components
    # --> the weights are applied on the Gram matrices
    
    for(i in 1:PC){
      KU[,,i] <- ratios[i]*do.call(get(paste(kernelY$fam, "_hsic", sep="")), 
                                   list(x=U[,i], param=param[i]))
    }
    
  }else{
    
    # computation of Gram matrices of the principal components
    # --> no weights
    
    for(i in 1:PC){
      KU[,,i] <- do.call(get(paste(kernelY$fam, "_hsic", sep="")), 
                         list(x=U[,i], param=param[i]))
    }
    
  }
  
  # computation of the final Gram matrix
  KY <- 1
  if(kernelY$combi=="sum"){
    for(i in 1:PC){
      KY <- KY+KU[,,i]
    }
  }else{
    for(i in 1:PC){
      KY <- KY*KU[,,i]
    }
  }
  
  # final update
  kernelY$fam <- rep(kernelY$fam, PC)
  if(!(kernelY$position=="nowhere")) kernelY$ratios <- ratios
  paramY <- param

  res <- list(KY=KY, kernelY=kernelY, paramY=paramY)
  
  return(res)
  
}

### GAK.KY #####################################################################

GAK.KY <- function(Y, paramY){
  
  ### inputs ###
  
  # Y:          n x q     matrix        (output samples)
  # paramY:     1 x 2     numeric       (parameter values for the GA kernel)
  
  ### output ###
  
  # KY:         n x n     matrix        (output Gram matrix)
  
  KY <- 1-as.matrix(proxy::dist(Y, method="gak", 
                                sigma=paramY[1], window.size=paramY[2], 
                                diag=TRUE, upper=TRUE))
  
  return(KY)
  
}

### DTW.KY #####################################################################

DTW.KY <- function(Y, kernelY){
  
  ### inputs ###
  
  # Y:          n x q     matrix        (output samples)
  # kernelY:    1 x ...   list          (options describing the kernel)
  
  ### output ###
  
  # KY:         n x n     matrix        (output Gram matrix)
  
  n <- nrow(Y)
  q <- ncol(Y)
  
  # a row-by-row split is operated on Y
  # --> the resulting object is a list of numeric (corresponding to row vectors)
  lot <- split(Y, rep(1:n, q))
  
  # the function "dtw_dismat" (from the IncDTW package) allows for much faster computations
  dist.dtw <- as.matrix(IncDTW::dtw_dismat(lot=lot, dist_metho="norm1")$dismat)
  
  # normalization
  nz.dtw <- dist.dtw/mean(dist.dtw)
  
  # the selected 1D-kernel is finally applied
  KY <- do.call(get(paste(kernelY$fam, "_hsic", sep="")),
                list(x=NULL, param=1, d=nz.dtw))
  
  return(KY)
  
}

### available.material #########################################################

available.material <- function(){
  
  ### no input ###

  ### output ###
  
  # L is a list containing the following elements:
  
  # kernels:            1 x ...   character   (available kernels)
  # anova.kernels:      1 x ...   character   (available anova kernels)
  # sd.kernels:         1 x ...   character   (kernels with one parameter of type sd)
  # free.kernels:       1 x ...   character   (parameter-free kernels)
  # weight.functions:   1 x ...   character   (available weight functions)
  # tests               1 x ...   character   (available independence tests)
  
  L <- list()
  
  ###############
  ### Kernels ###
  ###############
  
  ### Are they ANOVA kernels? ###
  
  # kernels that can be selected by the user in kernelX and kernelY
  
  L$kernels <- c("categ", "dcov", "invmultiquad",
                 "laplace", "laplace_anova", "linear", "matern3", "matern3_anova",
                 "matern5", "matern5_anova", "raquad", "rbf", "rbf_anova",
                 "sobolev1", "sobolev2")
  
  # kernels that directly allow for the HSIC-ANOVA decomposition
  
  L$anova.kernels <- c("sobolev1", "sobolev2",
                       "laplace_anova", "matern3_anova", "matern5_anova", "rbf_anova")
  
  # kernels that can be converted into ANOVA kernels
  
  L$conv.kernels <- c("laplace", "matern3", "matern5", "rbf")
  
  # kernels that can be coupled with PCA
  
  L$pca.kernels <- c("dcov", "invmultiquad", "laplace", "linear", "matern3", 
                     "matern5", "raquad", "rbf")
  
  # kernels that can be coupled with DTW
  
  L$dtw.kernels <- c("invmultiquad", "laplace", "matern3", 
                     "matern5", "raquad", "rbf")
  
  ### How to choose kernel parameters? ###

  # one-parameter kernels

  L$sd.kernels <- c("invmultiquad", "laplace", 
                    "matern3", "matern5", "raquad", "rbf", 
                    "laplace_anova", "matern3_anova", "matern5_anova", "rbf_anova")
  
  # parameter-free kernels
  
  L$free.kernels <- c("linear", "sobolev1", "sobolev2")
  
  ########################
  ### Weight functions ###
  ########################
  
  L$weight.functions <- c("indicTh", "zeroTh", "logistic", "exp1side")
  
  #############################
  ### Tests of independence ###
  #############################
  
  # test procedures that can be selected by the user in test.method
  
  L$tests <- c("Asymptotic", "Permutation", "Seq_Permutation", "Gamma")
  
  # test procedures that can be coupled with U-statistics
  
  L$U.tests <- c("Asymptotic", "Permutation", "Seq_Permutation")
  
  # test procedures that can be coupled with C-HSIC indices
  
  L$cond.tests <- c("Asymptotic", "Permutation", "Seq_Permutation", "Gamma")
  
  return(L)
  
}

### Kernels ####################################################################

categ_hsic <- function(x, param){
  if(param==0){
    d <- 1*sapply(x, function(z){z==x})
  }else{
    n <- length(x)
    d <- matrix(0, n, n)
    val <- unique(x)
    nval <- length(val)
    for(i in 1:nval){
      id <- which(x==val[i])
      d[id,id] <- 1/sum(x==val[i])
    }
  }
  return(d)
}

dcov_hsic <- function(x, param){
  n <- length(x)
  d <- as.matrix(dist(x))
  nrm <- matrix(abs(x), n, n)
  return(0.5*(nrm^param+t(nrm)^param-d^param))
}

invmultiquad_hsic <- function(x, param, d=NULL){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  d <- d/param
  return(1/sqrt(d^2+1))
}

laplace_hsic <- function(x, param, d=NULL){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  d <- d/param
  return(exp(-d))
}

linear_hsic <- function(x, ...){
  return(x%*%t(x))
}

matern3_hsic <- function(x, param, d=NULL){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  d <- d/param
  return((1+sqrt(3)*d)*exp(-sqrt(3)*d))
}

matern5_hsic <- function(x, param, d=NULL){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  d <- d/param
  return((1+sqrt(5)*d+5/3*d^2)*exp(-sqrt(5)*d))
}

raquad_hsic <- function(x, param, d=NULL){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  d <- d/param
  return(1-d^2/(1+d^2))
}

rbf_hsic <- function(x, param, d=NULL){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  d <- d/param
  return(exp(-0.5*d^2))
}

sobolev1_hsic <- function(x, ...){
  d <- as.matrix(dist(x))
  return(1+B1(x)%*%t(B1(x))+0.5*B2(d))
}

sobolev2_hsic <- function(x, ...){
  d <- as.matrix(dist(x))
  return(1+B1(x)%*%t(B1(x))+B2(x)%*%t(B2(x))/4-B4(d)/24)
}

### ANOVA kernels ##############################################################

# Ginsbourger, Roustant, Schuhmacher, Durrande & Lenz (2014)
# On ANOVA decompositions of kernels and Gaussian random field paths
# --> see Section 9: Additional examples

### formulas for the Laplace kernel ###

laplace_int <- function(x, param){
  k0 <- laplace_hsic(x=NULL, param=param, d=abs(x))
  k1 <- laplace_hsic(x=NULL, param=param, d=abs(1-x))
  vec <- param*(2-k0-k1)
  return(vec)
}

laplace_iint <- function(param){
  scl <- 2*param*(1-param+param*exp(-1/param))
  return(scl)
}

### formulas for the Matern 3/2 kernel ###

matern3_int <- function(x, param){
  param <- param/sqrt(3)
  alpha <- x/param
  beta <- (1-x)/param
  vec <- param*(4-( (2+alpha)*exp(-alpha) + (2+beta)*exp(-beta) ))
  return(vec)
}

matern3_iint <- function(param){
  param <- param/sqrt(3)
  scl <- 2*param*( 2-3*param+(1+3*param)*exp(-1/param) )
  return(scl)
}

### formulas for the Matern 5/2 kernel ###

matern5_int <- function(x, param){
  param <- param/sqrt(5)
  alpha <- x/param
  beta <- (1-x)/param
  vec <- (param/3)*(16 - (8+5*alpha+alpha^2)*exp(-alpha) 
                    - (8+5*beta+beta^2)*exp(-beta) )
  return(vec)
}

matern5_iint <- function(param){
  param <- param/sqrt(5)
  scl <- (param/3)*(16-30*param)+(2/3)*(1+7*param+15*param^2)*exp(-1/param)
  return(scl)
}

### formulas for the Gaussian kernel ###

rbf_int <- function(x, param){
  vec <- param*sqrt(2*pi)*(pnorm((1-x)/param)+pnorm(x/param)-1)
  return(vec)
}

rbf_iint <- function(param){
  scl <- 2*(param^2)*(exp(-1/(2*param^2))-1)+param*sqrt(2*pi)*(2*pnorm(1/param)-1)
  return(scl)
}

### kernel transformation ###

anova_transfo <- function(x, param, ker){
  
  ### inputs ###
  
  # x:        1 x n       numeric     (samples)
  # param:    1 x 1       numeric     (kernel parameter)
  # ker:      1 x 1       character   (kernel family)
  
  ### output ###
  
  # K0:       n x n       numeric     (Gram matrix)

  n <- length(x)
  
  # functions depending on ker
  ker.fun <- get(paste0(ker, "_hsic"))
  integ1.fun <- get(paste0(ker, "_int"))
  integ2.fun <- get(paste0(ker, "_iint"))
  
  #  initial Gram matrix
  K <- do.call(ker.fun, list(x=x, param=param))
  
  # formulas for integral calculus
  iK <- matrix(do.call(integ1.fun, list(x=x, param=param)), n, n)
  iiK <- do.call(integ2.fun, list(param=param))
  
  # new Gram matrix
  K0 <- 1+K-iK-t(iK)+iiK
  
  return(K0)
  
}

laplace_anova_hsic <- function(x, param){ anova_transfo(x, param, "laplace") }

matern3_anova_hsic <- function(x, param){ anova_transfo(x, param, "matern3") }

matern5_anova_hsic <- function(x, param){ anova_transfo(x, param, "matern5") }

rbf_anova_hsic <- function(x, param){ anova_transfo(x, param, "rbf") }

### Bernoulli polynomials ######################################################

B1 <- function(t){
  return(t-.5)
}

B2 <- function(t){
  return(t^2-t+1/6)
}

B4 <- function(t){
  return(t^4-2*t^3+t^2-1/30)
}

### print.sensiHSIC ############################################################

print.sensiHSIC <- function(x, ...){
  
  ### input ###
  
  # x:      1 x 1     list of class sensiHSIC  
  
  ##########
  
  # initial call:
  cat("\n","Call:","\n", deparse(x$call), "\n\n", sep="")
  
  if(!is.null(x$X)){
    # nb of model runs
    cat("Model runs: ", nrow(x$X), "\n\n", sep="")
    # R2 HSIC indices
    if(!is.null(x$S)){
      cat("R2-HSIC indices:", "\n", sep="")
      print(x$S)
      cat("\n")
    }
    # HSIC indices
    if(!is.null(x$HSICXY)){
      cat("HSIC indices:", "\n", sep="")
      print(x$HSICXY)
      cat("\n")
    }
    # first-order HSIC-ANOVA indices 
    if(!is.null(x$FO)){
      cat("First-order HSIC-ANOVA indices:", "\n", sep="")
      print(x$FO)
      cat("\n")
    }
    # total-order HSIC-ANOVA indices
    if(!is.null(x$TO)){
      cat("Total-order HSIC-ANOVA indices:", "\n", sep="")
      print(x$TO)
      cat("\n")
    }
  }else{
    cat("(empty)", "\n\n")
  }
  
}

### plot.sensiHSIC #############################################################

plot.sensiHSIC <- function(x, ylim=c(0, 1), ...){
  
  ### inputs ###
  
  # x:      1 x 1     list of class sensiHSIC 
  # ylim:   1 x 2     (numeric)     bounds on the y-axis
  
  ##########
  
  nb.stats <- sum(c("S", "FO", "TO") %in% names(x))
  
  if(nb.stats==0){
    
    warning("No result is stored in x. \n", call.=TRUE)
    
  }else{
    
    # nb of input variables
    p <- nrow(x$S) 
    
    ###################
    ### Preparation ###
    ###################
    
    xshift <- 0.15
    yshift <- 0.2
    
    # minimum and maximum values on the y-axis
    ym <- ylim[1]-yshift
    yM <- ylim[2]+nb.stats*yshift
    
    # points in the x-axis (to avoid overplot)
    xstart <- xval <- 1:p
    xfinal <- (1:p)+nb.stats*xshift
    
    # minimum and maximum values on the x-axis
    xm <- xstart[1]
    xM <- xfinal[p]
    
    # preliminary work for the legend
    
    cap.S <- "R2-HSIC indices"
    cap.FO <- "First-order HSIC-ANOVA indices"
    cap.TO <- "Total-order HSIC-ANOVA indices"
    
    marker.S <- 21
    marker.FO <- 22
    marker.TO <- 24
    
    col.S <- "blue"
    col.FO <- "red"
    col.TO <- "green"
    
    all.captions <- all.markers <- all.colors <- c()
    
    #######################
    ### R2-HSIC indices ###
    #######################
    
    if(!is.null(x$S)){
      
      nodeplot(x$S, at=xval,
               xlim=c(xm, xM), ylim=c(ym, yM),
               pch=marker.S, bg=col.S)
      
      xval <- xval+xshift
      
      all.captions <- c(all.captions, cap.S)
      all.markers <- c(all.markers, marker.S)
      all.colors <- c(all.colors, col.S)
      
    }
    
    #####################################################
    ### First-order HSIC-ANOVA indices (if available) ###
    #####################################################
    
    if(!is.null(x$FO)){
      
      nodeplot(x$FO, at=xval,
               xlim=c(xm, xM), ylim=c(ym, yM),
               pch=marker.FO, bg=col.FO,
               add=TRUE, labels=FALSE)
      
      xval <- xval+xshift
      
      all.captions <- c(all.captions, cap.FO)
      all.markers <- c(all.markers, marker.FO)
      all.colors <- c(all.colors, col.FO)
      
    }
    
    #####################################################
    ### Total-order HSIC-ANOVA indices (if available) ###
    #####################################################
    
    if(!is.null(x$TO)){
      
      nodeplot(x$TO, at=xval,
               xlim=c(xm, xM), ylim=c(ym, yM),
               pch=marker.TO, bg=col.TO,
               add=TRUE, labels=FALSE)
      
      xval <- xval+xshift
      
      all.captions <- c(all.captions, cap.TO)
      all.markers <- c(all.markers, marker.TO)
      all.colors <- c(all.colors, col.TO)
      
    }
    
    ##############
    ### Legend ###
    ##############
    
    abline(h=0, col="black", lwd=2, lty=2)
    abline(h=1, col="black", lwd=2, lty=2)
    
    legend("top", legend=all.captions, col="black", pch=all.markers, pt.bg=all.colors)
    
  }
  
}