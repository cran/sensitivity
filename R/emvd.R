############################################
# Needed packages
#library(parallel)
#library(gtools)
#library(boot)

############################################
#R2 estimation
estim.R2<-function(dataX,y, subset, logistic, any.cat){
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

    #Subsetting X
    X_<-dataX[,as.vector(subset)]
    dat<-data.frame(y,X_)
    sum_glm<-summary(glm(y~., data=dat, family="binomial"))

    #Computing Logistic R2 estimate
    R2<-1-(sum_glm$deviance/sum_glm$null.deviance)
  }
  return(R2)
}

######################################
#Recursive E-MVD computation
calc.emvd.rec.boot<-function(X, y, d, indices, logistic, any.cat, i=1:nrow(X)){
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

  #Bootstrap row selection
  X<-X[i,]
  y<-y[i,]

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
                        any.cat=any.cat)
  }

  #################################
  # E-MVD indices computation

  #Finding inputs with w({i})=0 (i.e., R2(full model) - R2(full model\i)=0)
  idx_null<-which(R2s[[d+1]]-R2s[[d]]==0)

  if(length(idx_null)!=0){
    #######################################
    #If there are any inputs with w({i})=0

    #Number of inputs with w({i})=0
    nv_null<-length(idx_null)
    #Position of inputs with w({i})=0
    v_null<-rev(1:d)[idx_null]
    #New number of dimension
    d_n<-d-nv_null
    #Preparing the recursive computation result, ommiting the excluded covariates
    Ps<-rep(list(0), d_n+1)
    #Same thing for the indices reference object
    indices_<-rep(list(0), d_n+1)
    for(j in 1:d_n){
      indices_[[j+1]]<-t(gtools::combinations(n=d_n, r=j))
      Ps[[j+1]]<-matrix(NA, ncol=choose(n=d_n, k=j), nrow=1)
    }
  }else{
    ###############################
    #If no excluded covariates
    indices_<-indices
    Ps<-rep(list(0), d+1)
    for(j in 1:d){
      Ps[[j+1]]<-matrix(NA, ncol=choose(n=d, k=j), nrow=1)
    }
    d_n<-d
  }

  #Preparing the recursive computation of the E-MVD indices
  #Full model indices
  N=1:d_n
  #Setting up w({emptyset})=1
  Ps[[1]]=1
  #Value of R2(full_model) to avoid a lot of search
  R2_comp<-R2s[[d+1]]

  for (ord in 1:(length(Ps)-1)){
    #For every input orders
    for(j in 1:ncol(Ps[[ord+1]])){
      #For every element of Ps, in each order
      #Find the corresponding subset
      S<-c(indices_[[ord+1]][,j])
      #Find the indices of the subset full_model\S
      NpS<-setdiff(N,S)
      if(ord==1){
        #If the subset is of size 1, then find the index of NpS in indices
        idx_nps<-which(apply(indices_[[d_n-ord+1]], 2,
                             function(x) all(x %in% c(NpS))))
        #Its P(.) value is initialized to R2(Full model) - R2(NpS)
        Ps[[ord + 1]][,j]=R2_comp - R2s[[d_n-ord+1]][idx_nps]
      }else if(ord==d_n){
        #If the subset is of size d, then w(S)=R2(full model)
        #Computing its P(.) value according to Feldman (2005)
        Ps[[ord+1]]=R2_comp*(1/sum(1/Ps[[ord]]))
      }else{
        #If the size of the subset is between 0 and d
        #Find the index of NpS in indices
        idx_nps<-which(apply(indices_[[d_n-ord+1]], 2,
                             function(x) all(x %in% c(NpS))))
        #Matrix of the S\i (on each column) for all i in S
        Spi<-matrix(sapply(S, function(x) setdiff(S, x)), ncol=length(S))
        #Vector of the positions of all the Spi in indices
        idx_Spi<-rep(0, ncol(Spi))
        for(k in 1:ncol(Spi)){
          idx_Spi[k]<-which(apply(indices_[[ord]],2,
                                  function(x) all(x %in%c(Spi[,k]))))
        }
        #Value of w(S)= R2(full model) - R2(full model\S)
        wS=R2_comp - R2s[[d-ord+1]][idx_nps]
        #Compute P(S) according to Feldman(2005)
        Ps[[ord+1]][,j]=wS*(1/sum(1/Ps[[ord]][,idx_Spi]))
      }
    }
  }
  #Prepare the output of the function
  if(length(idx_null)==0){
    #If there aren't any ommited covariate
    res<-matrix(rev(Ps[[d+1]]/Ps[[d]]), nrow=1)
  }else{
    #If there is any, setting up their E-MVD index to 0
    res_<-rep(0,d)
    res_[setdiff(1:d, v_null)]<-rev(Ps[[d_n+1]]/Ps[[d_n]])
    res<-matrix(res_, nrow=1)
  }
  colnames(res)=colnames(X)
  return(res)
}


############################################
# Main function
emvd<-function(X, y, logistic = FALSE,  rank=FALSE, nboot = 0, conf=0.95, parl=NULL){
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
                                       any.cat=any.cat)
      }else{
        #############################
        #R2 estimation
        R2s[[j+1]]<-apply(indices[[j+1]], 2,estim.R2,
                          dataX=X,
                          y=y,
                          logistic=logistic,
                          any.cat=any.cat)
      }
    }

    #################################
    # E-MVD indices computation

    #Finding inputs with w({i})=0 (i.e., R2(full model) - R2(full model\i)=0)
    idx_null<-which(R2s[[d+1]]-R2s[[d]]==0)

    if(length(idx_null)!=0){
      #######################################
      #If there are any inputs with w({i})=0

      #Number of inputs with w({i})=0
      nv_null<-length(idx_null)
      #Position of inputs with w({i})=0
      v_null<-rev(1:d)[idx_null]
      #New number of dimension
      d_n<-d-nv_null
      #Preparing the recursive computation result, ommiting the excluded covariates
      Ps<-rep(list(0), d_n+1)
      #Same thing for the indices reference object
      indices_<-rep(list(0), d_n+1)
      for(j in 1:d_n){
        indices_[[j+1]]<-t(gtools::combinations(n=d_n, r=j))
        Ps[[j+1]]<-matrix(NA, ncol=choose(n=d_n, k=j), nrow=1)
      }
    }else{
      ###############################
      #If no excluded covariates
      indices_<-indices
      Ps<-rep(list(0), d+1)
      for(j in 1:d){
        Ps[[j+1]]<-matrix(NA, ncol=choose(n=d, k=j), nrow=1)
      }
      d_n<-d
    }


    #Preparing the recursive computation of the E-MVD indices
    #Full model indices
    N=1:d_n
    #Setting up w({emptyset})=1
    Ps[[1]]=1
    #Value of R2(full_model) to avoid searching for it a lot
    R2_comp<-R2s[[d+1]]

    for (ord in 1:(length(Ps)-1)){
      #For every input orders
      for(j in 1:ncol(Ps[[ord+1]])){
        #For every element of Ps, in each order
        #Find the corresponding subset
        S<-c(indices_[[ord+1]][,j])
        #Find the indices of the subset full_model\S
        NpS<-setdiff(N,S)
        if(ord==1){
          #If the subset is of size 1, then find the index of NpS in indices
          idx_nps<-which(apply(indices_[[d_n-ord+1]], 2,
                               function(x) all(x %in% c(NpS))))
          #Its P(.) value is initialized to R2(Full model) - R2(NpS)
          Ps[[ord + 1]][,j]=R2_comp - R2s[[d_n-ord+1]][idx_nps]
        }else if(ord==d_n){
          #If the subset is of size d, then w(S)=R2(full model)
          #Computing its P(.) value according to Feldman (2005)
          Ps[[ord+1]]=R2_comp*(1/sum(1/Ps[[ord]]))
        }else{
          #If the size of the subset is between 0 and d
          #Find the index of NpS in indices
          idx_nps<-which(apply(indices_[[d_n-ord+1]], 2,
                               function(x) all(x %in% c(NpS))))
          #Matrix of the S\i (on each column) for all i in S
          Spi<-matrix(sapply(S, function(x) setdiff(S, x)), ncol=length(S))
          #Vector of the positions of all the Spi in indices
          idx_Spi<-rep(0, ncol(Spi))
          for(k in 1:ncol(Spi)){
            idx_Spi[k]<-which(apply(indices_[[ord]],2,
                                    function(x) all(x %in%c(Spi[,k]))))
          }
          #Value of w(S)= R2(full model) - R2(full model\S)
          wS=R2_comp - R2s[[d-ord+1]][idx_nps]
          #Compute P(S) according to Feldman(2005)
          Ps[[ord+1]][,j]=wS*(1/sum(1/Ps[[ord]][,idx_Spi]))
        }
      }
    }

    if(length(idx_null)==0){
      #If there aren't any ommited covariates
      emvd<-matrix(rev(Ps[[d+1]]/Ps[[d]]), ncol=1)
      rownames(emvd)=colnames(X)
      colnames(emvd)=c("original")
    }else{
      #If there is any, setting up their E-MVD index to 0
      res<-rep(0,d)
      res[setdiff(1:d, v_null)]<-rev(Ps[[d_n+1]]/Ps[[d_n]])
      emvd<-matrix(res, ncol=1)
      rownames(emvd)=colnames(X)
      colnames(emvd)=c("original")
    }

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
                                       any.cat=any.cat)
      }else{
        #############################
        #R2 estimation
        R2s[[j+1]]<-apply(indices[[j+1]], 2,estim.R2,
                          dataX=X,
                          y=y,
                          logistic=logistic,
                          any.cat=any.cat)
      }
    }

    if(is.null(parl)){
      #Bootstrap confidence interval estimation
      boot_emvd<-boot::boot(data=cbind(X),
                      statistic=calc.emvd.rec.boot,
                      R=nboot,
                      d=d,
                      y=y,
                      logistic=logistic,
                      indices=indices,
                      any.cat=any.cat,
                      stype="i")
    }else{
      #Parallelized Bootstrap confidence interval estimation
      boot_emvd<-boot::boot(data=cbind(X),
                     statistic=calc.emvd.rec.boot,
                     R=nboot,
                     d=d,
                     y=y,
                     logistic=logistic,
                     indices=indices,
                     any.cat=any.cat,
                     stype="i",
                     parallel="snow",
                     ncpus=parl,
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
