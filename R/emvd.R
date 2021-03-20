############################################
# Needed packages
#library(parallel)
#library(gtools)
#library(boot)

############################################
#R2 estimation
estim.R2<-function(dataX,y, subset, logistic, any.cat, max.iter){
  # Function estim.R2:
  #   Computes the R2 (explained variance) statistic for linear/logistic
  #   regression for a specified sub-model (subset of inputs)
  # Inputs:
  #   dataX: data.frame/matrix
  #   y: num/factor/logical vector
  #   subset: integer vector, subset of inputs for which to compute the R2
  #   logistic: logical. If TRUE, then logistic regression is performed
  #   any.cat: logical. If TRUE, then at least one column of dataX is a factor.
  # Output:
  #   R2: num. R2 value for the desired model (linear/logistic) and for the
  #       specified subset of inputs.

  if(!logistic){
    #######################
    #Linear Model

    if(!any.cat){
      ################################################################
      # Computing by hand with scalar inputs: approx 7x faster than lm

      #Subsetting X
      X_<-as.matrix(dataX[,as.vector(subset)])

      #Adding beta_0
      X_<-cbind(rep(1, nrow(X_)), X_)

      #Computing beta estimates
      estim.beta<-solve(t(X_)%*%X_)%*%t(X_)%*%y

      #Computing R2 estimate
      R2<-1-(var(y-X_%*%estim.beta)/var(y))
    }else{
      #################################################################
      # Using lm method if there is any categorical inputs

      #Subsetting X
      X_<-dataX[,as.vector(subset)]
      dat<-data.frame(y,X_)

      #Computing R2 estimate
      R2<-summary(lm(y~., data=dat))$r.squared
    }
  }else{
    #######################
    #Logistic Model
    #Maximum optimization iterations
    fit.control=glm.control(maxit=max.iter)
    #Subsetting X
    X_<-dataX[,as.vector(subset)]
    dat<-data.frame(y,X_)
    sum_glm<-summary(glm(y~., data=dat, 
                         family="binomial",
                         control=fit.control))

    #Computing Logistic R2 estimate
    R2<-1-(sum_glm$deviance/sum_glm$null.deviance)
  }
  return(R2)
}

######################################
#Recursive E-MVD computation
calc.emvd.rec.boot<-function(X, y, indices, logistic, max.iter,tol, any.cat, i=1:nrow(X)){
  # Function calc.emvd.rec.boot
  #   Computes the E-MVD indices using the recursive algorithm of Feldman (2005)
  #   for bootstrap confidence interval computations.
  # Inputs:
  #   X: matrix of data frame.
  #   y: num/factor/logical vector
  #   d: int. Covariates dimension.
  #   indices: list. List of subset indices.
  #   logistic: logical. If TRUE, logistic regression is performed, else, linear
  #             model instead.
  #   any.cat: logical. If TRUE, then at least one column of dataX is a factor.
  #   i: integer vector. Used for bootstrap estimation.
  # Output:
  #   res: (1xd) matrix. E-MVD value for all covariates.

  #################################
  # Pre-Processing
  #Dimension
  d=ncol(X)
  #Bootstrap row selection
  X<-X[i,]
  y<-matrix(y,ncol=1)[i,]

  #Initiating the conditional elements object
  R2s<-rep(list(0), d+1)

  #################################
  # Conditional elements estimation
  for(j in 1:d){
      #############################
      #R2 estimation
      R2s[[j+1]]<-apply(indices[[j+1]], 2,estim.R2,
                        dataX=X,
                        y=y,
                        logistic=logistic,
                        any.cat=any.cat,
                        max.iter=max.iter)
  }

  #################################
  # E-MVD indices computation
  
  #Preparing the recursive computation of the E-MVD indices
  #Full model indices
  N=1:d
  #Initiating P(S)
  Ps<-rep(list(0), d+1)
  #Computing w({i}) forall i
  Ps[[1]]=rev(R2s[[d+1]] - R2s[[d]])
  
  #Finding if any w({i})=0 or |w({i})|<tol
  if(is.null(tol)){
    logi.zeros<-Ps[[1]]==0
  }else{
    logi.zeros<-abs(Ps[[1]])<tol
  }
  
  if(any(logi.zeros)){
    var_null<-which(logi.zeros)
    nb_null<-length(var_null)
    idx_ind_null=rep(list(0), d-nb_null+1)
    indices_<-rep(list(0), d-nb_null)
    for(i in 1:(d-nb_null)){
      checkmat<-matrix(is.element(indices[[i+1]], var_null), i, ncol(indices[[i+1]]))
      idx_ind_null[[i]] <- as.vector(which(colSums(checkmat) == 0))
      indices_[[i]]<-matrix(indices[[i+1]][, idx_ind_null[[i]]],ncol=length(idx_ind_null[[i]]))
    }
    p=d #Original number of dimensions
    d=d-nb_null #Number of dimensions after removing spurious covariates
    Ps[[1]]<-Ps[[1]][setdiff(N, var_null)] #Only keeping w({i}) for non-spurious covariates
  }else{
    p=d
    nb_null=0
  }
  #Value of R2(full_model) to avoid a lot of search
  R2_comp<-R2s[[p+1]]
  
  for(ord in 2:d){
    #For every input orders : we start at 2 because Ps[[1]] has already been computed
    if(ord==d){
      #If this is the last order
      Ps[[ord]]=R2_comp*(1/sum(1/Ps[[ord-1]]))
    }else{
      if(any(logi.zeros)){
        #########################
        #If any spurious variables
        #Computing w(S) for all S of order ord
        Ws<-rev(R2_comp - R2s[[d+1-ord]][idx_ind_null[[d-ord]]])
        nbset<-ncol(indices_[[ord]]) #Number of subsets S
        Ps[[ord]]<-rep(0, nbset) #Setting up the Ps results
        for(i in 1:nbset){
          S<-c(indices_[[ord]][,i]) #For every possible subset S of order ord
          idx_Spi=which(colSums(matrix(indices_[[ord-1]]%in%S, nrow=length(S)-1, ncol=ncol(indices_[[ord-1]])))==length(S)-1) #Find every S\{i} for all i
          Ps[[ord]][i]=Ws[i]/sum(1/(Ps[[ord-1]][idx_Spi])) #Recursively compute Ps
        }
        
      }else{
        ########################
        #If no spurious variables
        #Computing w(S) for all S of order ord
        Ws<-rev(R2_comp - R2s[[d+1-ord]])
        nbset<-ncol(indices[[ord+1]])#Number of subsets S
        Ps[[ord]]<-rep(0, nbset)#Setting up the Ps results
        for(i in 1:nbset){
          S<-c(indices[[ord+1]][,i]) #For every possible subset S of order ord
          idx_Spi=which(colSums(matrix(indices[[ord]]%in%S, nrow=length(S)-1, ncol=ncol(indices[[ord]])))==length(S)-1)#Find every S\{i} for all i
          Ps[[ord]][i]=Ws[i]/sum(1/(Ps[[ord-1]][idx_Spi]))#Recursively compute Ps
        }
      }
    }
  }
  
  pmd<-matrix(rev(Ps[[d]]/Ps[[d-1]]), nrow=1)
  if(any(logi.zeros)){
    #If spurious variable, we set their value to 0
    index <- 1
    pmd_ <- rep(0, p)
    for (i in 1:p) {
      if (!is.element(i, var_null)) {
        pmd_[i] <- pmd[index]
        index <- index + 1
      }
    }
    pmd=matrix(pmd_, nrow=1)
  }
  colnames(pmd)=colnames(X)
  return(pmd)
}


############################################
# Main function
emvd<-function(X, y, logistic = FALSE,  tol=NULL, rank=FALSE, nboot = 0, conf=0.95,  max.iter=1000, parl=NULL){
  ###########################################
  #Pre-processing

  #Retrieving dimension
  d<-ncol(X)

  #Checking for categorical inputs
  any.cat<-any(sapply(X,class)=="factor")

  ##########################################
  #Arguments check

  #X
  if(!(class(X)[1]=="matrix" | class(X)[1]=="data.frame")){
    stop("X must be either a matrix of a data frame.")
  }

  #y
  if (logistic==F & !is.numeric(y)){
    stop(paste("y must be a numeric vector of length ",nrow(X)," or a matrix with ", nrow(X), " rows and one column.", sep=""))
  }else if(logistic==T & (all(!is.logical(y), !is.factor(y), !is.numeric(y))) ){
    stop("y must be a logical, factor or numeric vector.")
  }else if (length(y)!=nrow(X)){
    stop(paste("y must be a vector of length ",nrow(X)," or a matrix with ", nrow(X), " rows and one column.", sep=""))
  }

  #logistic
  if(!is.logical(logistic)){
    stop("The 'logistic' argument must be logical, either TRUE or FALSE.")
  }

  #rank
  if(!is.logical(rank)){
    stop("The 'rank' argument must be logical, either TRUE or FALSE.")
  }else if (logistic & rank){
    rank=F
    warning("Impossible to perform a logistic regression with a rank transformation. Defaulted to rank=FALSE.")
  }else if (any.cat & rank){
    rank=F
    warning("Impossible to perform a rank transformation with categorical inputs. Defaulted to rank=FALSE.")
  }

  #nboot
  if(!is.numeric(nboot)){
    stop("The 'nboot' argument must be a positive integer.")
  }else if (nboot%%1!=0 | nboot<0 | length(nboot)!=1){
    stop("The 'nboot' argument must be a positive integer.")
  }else if(nboot==0){
    boot=F
  }else{
    boot=T
  }

  #parl
  if(!is.numeric(parl)){
    if(!is.null(parl)){
      stop("The 'parl' argument must be NULL or a positive integer.")
    }
  }else if(parl%%1!=0 | parl<0 | length(parl)!=1){
    stop("The 'parl' argument must be NULL or a positive integer.")
  }else if (parl>parallel::detectCores()){
    parl=parallel::detectCores()-1
    warning(paste("Too many cores specified. Defaulted to ", parl, " cores.", sep=""))
  }

  ########################################
  #Data Processing

  #Rank Transformation
  if(rank){
    X<-apply(X, 2, rank)
  }

  #List of subsets by order ([[1]]: subsets of order 1...etc)
  indices<-rep(list(0), d+1)
  #List of R2s, arranged by orders (same structure as the indices object)
  R2s<-rep(list(0), d+1)
  #List of elements for recursive computing

  #Computing subsets
  for(j in 1:d){
    indices[[j+1]]<-t(gtools::combinations(n=d, r=j))
  }

  #####################################################
  # E-MVD Computing

  #Preparing parallel cluster
  if(!is.null(parl)){
    cl=parallel::makeCluster(parl)
  }

  if(!boot){
    #################################
    # Non-Bootstrap computation

    ############################################
    #Conditional elements computing
    #For every order, the corresponding R2 estimate of each
    #subset of indices are stored in R2s
    for(j in 1:d){
      if(!is.null(parl)){
        #############################
        #Parallelized R2 estimation
        R2s[[j+1]]<-parallel::parApply(cl,indices[[j+1]], 2,estim.R2,
                                       dataX=X,
                                       y=y,
                                       logistic=logistic,
                                       any.cat=any.cat,
                                       max.iter=max.iter)
      }else{
        #############################
        #R2 estimation
        R2s[[j+1]]<-apply(indices[[j+1]], 2,estim.R2,
                          dataX=X,
                          y=y,
                          logistic=logistic,
                          any.cat=any.cat,
                          max.iter=max.iter)
      }
    }

    #################################
    # E-MVD indices computation
    #Preparing the recursive computation of the E-MVD indices
    #Full model indices
    N=1:d
    #Initiating P(S)
    Ps<-rep(list(0), d+1)
    #Computing w({i}) forall i
    Ps[[1]]=rev(R2s[[d+1]] - R2s[[d]])
    
    #Finding if any w({i})=0 or |w({i})|<tol
    if(is.null(tol)){
      logi.zeros<-Ps[[1]]==0
    }else{
      logi.zeros<-abs(Ps[[1]])<tol
    }
    
    if(any(logi.zeros)){
      var_null<-which(logi.zeros)
      nb_null<-length(var_null)
      idx_ind_null=rep(list(0), d-nb_null+1)
      indices_<-rep(list(0), d-nb_null)
      for(i in 1:(d-nb_null)){
        checkmat<-matrix(is.element(indices[[i+1]], var_null), i, ncol(indices[[i+1]]))
        idx_ind_null[[i]] <- as.vector(which(colSums(checkmat) == 0))
        indices_[[i]]<-matrix(indices[[i+1]][, idx_ind_null[[i]]],ncol=length(idx_ind_null[[i]]))
      }
      p=d #Original number of dimensions
      d=d-nb_null #Number of dimensions after removing spurious covariates
      Ps[[1]]<-Ps[[1]][setdiff(N, var_null)] #Only keeping w({i}) for non-spurious covariates
    }else{
      p=d
      nb_null=0
    }
    #Value of R2(full_model) to avoid a lot of search
    R2_comp<-R2s[[p+1]]
    
    if(!is.null(parl)){doParallel::registerDoParallel(cl)}
    for(ord in 2:d){
      #For every input orders : we start at 2 because Ps[[1]] has already been computed
      if(ord==d){
        #If this is the last order
        Ps[[ord]]=R2_comp*(1/sum(1/Ps[[ord-1]]))
      }else{
        if(any(logi.zeros)){
          #########################
          #If any spurious variables
          #Computing w(S) for all S of order ord
          Ws<-rev(R2_comp - R2s[[d+1-ord]][idx_ind_null[[d-ord]]])
          nbset<-ncol(indices_[[ord]]) #Number of subsets S
         
          if(!is.null(parl)){
            Ps[[ord]]<-foreach::foreach(i=1:nbset, .combine=cbind)%dopar%{
              S<-c(indices_[[ord]][,i])
              idx_Spi=which(colSums(matrix(indices_[[ord-1]]%in%S, nrow=length(S)-1, ncol=ncol(indices_[[ord-1]])))==length(S)-1)
              res=Ws[i]/sum(1/(Ps[[ord-1]][idx_Spi]))
              res
            }
          }else{
            Ps[[ord]]<-rep(0, nbset) #Setting up the Ps results
            for(i in 1:nbset){
              S<-c(indices_[[ord]][,i]) #For every possible subset S of order ord
              idx_Spi=which(colSums(matrix(indices_[[ord-1]]%in%S, nrow=length(S)-1, ncol=ncol(indices_[[ord-1]])))==length(S)-1) #Find every S\{i} for all i
              Ps[[ord]][i]=Ws[i]/sum(1/(Ps[[ord-1]][idx_Spi])) #Recursively compute Ps
            }
          }
        }else{
          ########################
          #If no spurious variables
          #Computing w(S) for all S of order ord
          Ws<-rev(R2_comp - R2s[[d+1-ord]])
          nbset<-ncol(indices[[ord+1]])#Number of subsets S
         
          if(!is.null(parl)){
            Ps[[ord]]<-foreach::foreach(i=1:nbset, .combine=cbind)%dopar%{
              S<-c(indices[[ord+1]][,i])
              idx_Spi=which(colSums(matrix(indices[[ord]]%in%S, nrow=length(S)-1, ncol=ncol(indices[[ord]])))==length(S)-1)
              res=Ws[i]/sum(1/(Ps[[ord-1]][idx_Spi]))
              res
            }
          }else{
            Ps[[ord]]<-rep(0, nbset)#Setting up the Ps results
            for(i in 1:nbset){
              S<-c(indices[[ord+1]][,i]) #For every possible subset S of order ord
              idx_Spi=which(colSums(matrix(indices[[ord]]%in%S, nrow=length(S)-1, ncol=ncol(indices[[ord]])))==length(S)-1)#Find every S\{i} for all i
              Ps[[ord]][i]=Ws[i]/sum(1/(Ps[[ord-1]][idx_Spi]))#Recursively compute Ps
            }
          }
        }
      }
    }
    
    emvd<-matrix(rev(Ps[[d]]/Ps[[d-1]]), ncol=1)
    if(any(logi.zeros)){
      #If spurious variable, we set their value to 0
      index <- 1
      emvd_ <- rep(0, p)
      for (i in 1:p) {
        if (!is.element(i, var_null)) {
          emvd_[i] <- emvd[index]
          index <- index + 1
        }
      }
      emvd=matrix(emvd_, ncol=1)
    }
    rownames(emvd)=colnames(X)
    colnames(emvd)=c("original")

    #Preparing the output.
    out<-list("call"=match.call(),
              "emvd"=emvd,
              "R2s"=R2s,
              "indices"=indices,
              "P"=Ps,
              "conf_int"=NULL,
              "X"=X,
              "y"=y,
              "logistic"=logistic,
              "boot"=boot,
              "nboot"=nboot,
              "rank"=rank,
              "parl"=parl,
              "conf"=conf)
  }else{
    ##################################
    # Bootstrap computation

    ############################################
    #Conditional elements computing
    #For every order, the corresponding R2 estimate of each
    #subset of indices are stored in R2s
    for(j in 1:d){
      if(!is.null(parl)){
        #############################
        #Parallelized R2 estimation
        R2s[[j+1]]<-parallel::parApply(cl,indices[[j+1]], 2,estim.R2,
                                       dataX=X,
                                       y=y,
                                       logistic=logistic,
                                       any.cat=any.cat,
                                       max.iter=max.iter)
      }else{
        #############################
        #R2 estimation
        R2s[[j+1]]<-apply(indices[[j+1]], 2,estim.R2,
                          dataX=X,
                          y=y,
                          logistic=logistic,
                          any.cat=any.cat,
                          max.iter=max.iter)
      }
    }

    if(is.null(parl)){
      #Bootstrap confidence interval estimation
      boot_emvd<-boot::boot(data=cbind(X),
                      statistic=calc.emvd.rec.boot,
                      R=nboot,
                      # d=d,
                      y=y,
                      logistic=logistic,
                      indices=indices,
                      any.cat=any.cat,
                      tol=tol,
                      max.iter=max.iter,
                      stype="i")
    }else{
      #Parallelized Bootstrap confidence interval estimation
      boot_emvd<-boot::boot(data=cbind(X),
                     statistic=calc.emvd.rec.boot,
                     R=nboot,
                     # d=d,
                     y=y,
                     logistic=logistic,
                     indices=indices,
                     any.cat=any.cat,
                     stype="i",
                     parallel="snow",
                     ncpus=parl,
                     tol=tol,
                     max.iter=max.iter,
                     cl=cl)

    }

    res.boot<-bootstats(boot_emvd, conf, "basic")
    rownames(res.boot)<-colnames(X)
    emvd<-matrix(res.boot[,1], ncol=1)
    colnames(emvd)=c("original")
    rownames(emvd)=colnames(X)
    out<-list("call"=match.call(),
              "emvd"=emvd,
              "R2s"=R2s,
              "indices"=indices,
              "P"=NULL,
              "conf_int"=res.boot,
              "X"=X,
              "y"=y,
              "logistic"=logistic,
              "boot"=boot,
              "nboot"=nboot,
              "rank"=rank,
              "parl"=parl,
              "conf"=conf)
  }

  #Stopping the cluster after parallelized computation
  if(!is.null(parl)){parallel::stopCluster(cl)}

  class(out)<-"emvd"

  return(out)
}

#######################################
# Custom methods
print.emvd<-function(x, ...){
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if(x$logistic){
    cat("\nE-MVD decomposition of R2 for logistic model\n")
  }else{
    cat("\nE-MVD decomposition of R2 for linear model\n")
  }
  if(x$boot){
    print(x$conf_int)
  }else{
    print(x$emvd)
  }
}


plot.emvd<-function(x, ylim=c(0,1), ...){
  if(x$logistic){
    plot_title="E-MVD decomposition of R2 for logistic model"
  }else{
    plot_title="E-MVD decomposition of R2 for linear model"
  }

  plot(x$emvd,
       ylim=ylim,
       axes=F,
       xlab="Inputs",
       ylab="E-MVD",
       main=plot_title,
       ...)
  axis(2)
  axis(1, at=seq_along(x$emvd), labels=colnames(x$X))
  box()
  graphics::grid()

  if(x$boot){
    segments(x0=1:ncol(x$X),
             y0=x$conf_int$`min. c.i.`,
             x1=1:ncol(x$X),
             y1=x$conf_int$`max. c.i.`)
    legend("topright",
           pch=c(1,NA),
           lty=c(NA,1),
           legend=c("E-MVD values",
                    paste(x$conf*100, "% Confidence interval", sep="")),
           bg="white")
  }else{
    legend("topright",
           pch=1,
           legend=c("E-MVD Values"),
           bg="white")
  }
}
