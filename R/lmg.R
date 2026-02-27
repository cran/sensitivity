############################################
# Needed packages
#library(parallel)
#library(doParallel)
#library(foreach)
#library(gtools)
#library(boot)

############################################
#R2 estimation
estim.R2<-function(dataX, y, subset, logistic, any.cat, max.iter){
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

############################################
# LMG computation function
calc.lmg<-function(X, y, d, logistic, indices, comb_weights, any.cat, parl,clus, boot,max.iter, i=1:nrow(X)){

  ###############################
  #Pre-Processing
  #Row selection (for Bootstrap estimate, default doesn't change anything)
  X<-X[i,]
  y<-y[i]

  #List of R2s, arranged by orders (same structure as the indices object)
  R2s<-rep(list(0), d+1)

  ####################################
  #Conditional elements estimation

  #For every order, the corresponding R2 estimate of each
  #subset of indices are stored in R2s
  for(j in 1:d){
    if(!is.null(parl)){
      #############################
      #Parallelized R2 estimation
      R2s[[j+1]]<-parallel::parApply(clus,indices[[j+1]], 2,estim.R2,
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
  ###############################
  #Affecting the weights
  if(!is.null(parl)){
    ################################
    #Parallelized weight affectation

    #One core per input
    res_lmg <- foreach::foreach(i=1:d, .combine=cbind)%dopar%{
      res=0
      #For every order, the increments are weighted and summed
      #for the input i
      for(ord in 1:d){
        if(ord==1){
          idx<-which(indices[[ord+1]]==i)
          res=res+comb_weights[ord]*R2s[[ord+1]][as.vector(idx)]
        }else{
          #Columns of subsets of order ord+1 containing var_j
          idx_j<-which(indices[[ord+1]]==i, arr.ind=T)[,2]
          #Columns of subsets of order ord without j
          idx_woj<-which(apply(indices[[ord]]!=i, 2, all))
          #Total increment value for order ord
          tot_incr<-sum(R2s[[ord+1]][as.vector(idx_j)]) - sum(R2s[[ord]][as.vector(idx_woj)])
          res=res+comb_weights[ord]*tot_incr
        }
      }
      res
    }
  }else{
    res_lmg<-rep(0,d)
    for(var_j in 1:d){
      #For every order, the increments are weighted and summed
      #for the input var_j
      for(ord in 1:d){
        if(ord==1){
          idx<-which(indices[[ord+1]]==var_j)
          res_lmg[var_j]=res_lmg[var_j]+comb_weights[ord]*R2s[[ord+1]][as.vector(idx)]
        }else{
          #Columns of subsets of order ord+1 containing var_j
          idx_j<-which(indices[[ord+1]]==var_j, arr.ind=T)[,2]
          #Columns of subsets of order ord without j
          idx_woj<-which(apply(indices[[ord]]!=var_j, 2, all))
          #Total increment value for order ord
          tot_incr<-sum(R2s[[ord+1]][as.vector(idx_j)]) - sum(R2s[[ord]][as.vector(idx_woj)])
          res_lmg[var_j]=res_lmg[var_j]+comb_weights[ord]*tot_incr
        }
      }
    }
  }

  #The results are divided by the dimension
  res_lmg=res_lmg/d

  if(boot){
    #If Bootstrap computation, return the LMG
    return(res_lmg)
  }else{
    #For initial computation, return conditional elements and LMG
    return(list("R2"=R2s,
                "lmg"=res_lmg))
  }
}

############################################
# Main function
lmg<-function(X, y, logistic = FALSE,  rank=FALSE, nboot = 0, conf=0.95, max.iter=1000, parl=NULL){
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
  #Weights vector, each element is the weight for |A| ([1]:|A|=0, [2]:|A|=1, ...)
  comb_weights<-rep(0,d)

  #Computing subsets and weights
  for(j in 1:d){
    indices[[j+1]]<-t(gtools::combinations(n=d, r=j))
    comb_weights[j]<-1/choose(n=(d-1), j-1)
  }

  ########################################
  #LMG Computation

  #Preparing parallel cluster
  if(!is.null(parl)){
    cl=parallel::makeCluster(parl)
    doParallel::registerDoParallel(cl)
  }else{
    cl=NULL
  }

  #Initial computation of the LMG indices
  res_lmg<-calc.lmg(X=X,
                   y=y,
                   d=d,
                   logistic=logistic,
                   indices=indices,
                   comb_weights=comb_weights,
                   any.cat=any.cat,
                   parl=parl,
                   clus=cl,
                   max.iter=max.iter,
                   boot=F)
  #Retrieving input names
  res_lmg$lmg<-matrix(res_lmg$lmg, ncol=1)
  colnames(res_lmg$lmg)<-"original"
  rownames(res_lmg$lmg)<-colnames(X)

  #Preparing the output
  out<-list("call"=match.call(),
            "lmg"=res_lmg$lmg,
            "R2s"=res_lmg$R2,
            "indices"=indices,
            "w"=comb_weights,
            "conf_int"=NULL,
            "X"=X,
            "y"=y,
            "logistic"=logistic,
            "boot"=boot,
            "nboot"=nboot,
            "rank"=rank,
            "parl"=parl,
            "conf"=conf
  )

  #Bootstrap C.I estimation
  if(boot){
    #Computing Boostrap C.Is
    if(is.null(parl)){
      boot_lmg<-boot(data=cbind(X),
                     statistic=calc.lmg,
                     R=nboot,
                     d=d,
                     y=y,
                     logistic=logistic,
                     indices=indices,
                     comb_weights=comb_weights,
                     any.cat=any.cat,
                     parl=NULL,
                     clus=NULL,
                     boot=boot,
                     max.iter=max.iter,
                     stype="i")
    }else{
      boot_lmg<-boot(data=cbind(X),
                     statistic=calc.lmg,
                     R=nboot,
                     d=d,
                     y=y,
                     logistic=logistic,
                     indices=indices,
                     comb_weights=comb_weights,
                     any.cat=any.cat,
                     parl=NULL,
                     clus=NULL,
                     boot=boot,
                     stype="i",
                     parallel="snow",
                     ncpus=parl,
                     max.iter=max.iter,
                     cl=cl)
    }

    #Retrieving input names
    CI_lmg<-bootstats(boot_lmg, conf, "basic")
    rownames(CI_lmg)=colnames(X)
    #Adding Bootstrap C.I estimations
    out$conf_int<-CI_lmg
  }
  #Stopping the cluster after parallelized computation
  if(!is.null(parl)){parallel::stopCluster(cl)}

  #Customize the class of the output for custom methods
  class(out)<-"lmg"

  return(out)
}


#######################################
# Custom methods

print.lmg<-function(x, ...){
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if(x$logistic){
    cat("\nDecomposition of R2 for logistic model\n")
  }else{
    cat("\nDecomposition of R2 for linear model\n")
  }
  if(x$boot){
    print(x$conf_int)
  }else{
    print(x$lmg)
  }
}


plot.lmg<-function(x, ylim=c(0,1), ...){
  if(x$logistic){
    plot_title="Decomposition of R2 for logistic model"
  }else{
    plot_title="Decomposition of R2 for linear model"
  }

  plot(x$lmg,
       ylim=ylim,
       axes=F,
       xlab="Inputs",
       ylab="LMG",
       main=plot_title,
       ...)
  axis(2)
  axis(1, at=seq_along(x$lmg), labels=colnames(x$X))
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
           legend=c("LMG values",
                    paste(x$conf*100, "% Confidence interval", sep="")),
           bg="white")
  }else{
    legend("topright",
           pch=1,
           legend=c("LMG Values"),
           bg="white")
    }
}
