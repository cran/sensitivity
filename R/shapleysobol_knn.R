############################################
# Needed packages
# library(parallel)
# library(doParallel)
# library(foreach)
# library(gtools)
# library(boot)
# library(RANN)
# library(whitening)
# library(TSP)

#########################
#Input rescaling function
rescale_inputs <- function(X){
  # Mahalanobis whitening, new X has unit diagonal covariance matrix
  X <- whitening::whiten(X, method="ZCA-cor")
  # We then apply a copula transform
  X <- apply(X,2,rank)/nrow(X)
  return(X)
}

############################################
# Var(E[Y|X_A])/Var(Y) estimation by Rank/TSP
estim.VE.rank<-function(dataX, y, subset, n){
  #Subsetting X
  X<-dataX[, subset, drop=F]
  any.cat<-any(sapply(X,class)=="factor")
  p=ncol(X)
  y=c(y)
  
  if(any.cat){
    id.cat=which(sapply(X,class)=="factor")
    id.cont <- setdiff(1:p,id.cat)
    Xtemp <- as.matrix(X[,id.cont,drop=FALSE])
    for (i in 1:length(id.cat)){
      X.factor<-X[,id.cat[i], drop=F]
      df.temp <- X[,id.cat[i],drop=FALSE]; colnames(df.temp) <- "X"
      # Normalization to ensure that two different rows will have unit euclidean distance
      xx.temp <- model.matrix(~X-1,df.temp)*sqrt(2)/2
      Xtemp<-cbind(Xtemp, xx.temp)
    }
    X<-Xtemp
  }
  
  # Use ranking or TSP
  if (ncol(X)>1){
    # Use TSP
    dist.mat <- as.matrix(dist(X))
    tsp <- TSP::TSP(dist.mat)
    tour <- as.integer(TSP::solve_TSP(tsp))
    id.sort <- cbind(tour,tour[c(2:n,1)])
  }else{
    # Use ranking
    id.sort <- sort(X[,1], index.return = TRUE)$ix
    id.sort <- matrix(c(id.sort,id.sort[c(2:n,1)]),ncol=2)
  }
  
  var.cond<-apply(matrix(y[id.sort],ncol=2),1,var)
  
  #Returning an estimate of Var(E[Y|X_A])/Var(Y)
  VE_A<-1-mean(var.cond)/var(y)
  return(VE_A)
}

############################################
# Var(E[Y|X_A])/Var(Y) estimation by KNN
estim.VE.knn<-function(dataX, y, subset, n, n.knn, n.limit){
  #Subsetting X
  X<-dataX[, subset, drop=F]
  any.cat<-any(sapply(X,class)=="factor")
  p=ncol(X)
  y=c(y)
  
  #One Hot Encoding of categorical variables
  if(any.cat){
    id.cat=which(sapply(X,class)=="factor")
    id.cont <- setdiff(1:p,id.cat)
    Xtemp <- as.matrix(X[,id.cont,drop=FALSE])
    for (i in 1:length(id.cat)){
      X.factor<-X[,id.cat[i], drop=F]
      df.temp <- X[,id.cat[i],drop=FALSE]; colnames(df.temp) <- "X"
      # Normalization to ensure that two different rows will have unit euclidean distance
      xx.temp <- model.matrix(~X-1,df.temp)*sqrt(2)/2
      Xtemp<-cbind(Xtemp, xx.temp)
    }
    X<-Xtemp
  }
  if(n<=n.limit){
    ###############################################
    # When n is reasonable, compute exact distances
    
    #Computing distance matrix
    dist.mat<-as.matrix(dist(X))
    #Getting the n.knn closest points for each datapoint
    id.sort <- t(apply(dist.mat,2,order)[1:n.knn,])
  }else{
    ###############################################
    # When n is too large, compute approximated NN
    
    #Getting the approximated n.knn closest points for each datapoint
    id.sort <- RANN::nn2(X,k=n.knn)$nn.idx
  }
  #Estimating Var(Y|X_A)
  var.cond<-apply(matrix(y[id.sort],ncol=n.knn),1,var)
  
  #Returning an estimate of Var(E[Y|X_A])/Var(Y)
  VE_A<-1-mean(var.cond)/var(y)
  return(VE_A)
}

###########################################
# Computing Shapley effects via KNN 
calc.shapknn<-function(X, y, method, n.perm, perms, d, indices, comb_weights, parl, clus, boot, n.limit, n.knn, noise, i=1:nrow(X), get.sob=F){
  ###############################
  #Pre-Processing
  #Row selection (for Bootstrap estimate, default doesn't change anything)
  if(boot){
    y<-X[,1]
    X<-X[,-1]
  }
  n<-nrow(X)
  
  #List of VEs, arranged by orders (same structure as the indices object)
  if(get.sob){
    VEs<-rep(list(0), length(indices))
    d<-length(indices)-1
  }else{
    VEs<-rep(list(0), d+1)
  }

  
  ####################################
  #Conditional elements estimation
  #For every order, the corresponding VE estimate of each
  #subset of indices are stored in VEs
  if(noise){
    #The Shapley sum up to V(E[Y|X_D])/V(Y) with D={1,...,d} (all variables)
    #which may not be equal to 1
    for(j in 1:d){
      if(method=='knn'){
        #############################################
        # KNN Computation of the conditional elements
        if(!is.null(parl)){
          #############################
          #Parallelized VE estimation
          VEs[[j+1]]<-parallel::parApply(clus,indices[[j+1]], 2,estim.VE.knn,
                                         dataX=X,
                                         y=y,
                                         n=n,
                                         n.knn=n.knn,
                                         n.limit=n.limit)
          
        }else{
          #############################
          # VE estimation
          VEs[[j+1]]<-apply(indices[[j+1]], 2, estim.VE.knn,
                            dataX=X,
                            y=y,
                            n=n,
                            n.knn=n.knn,
                            n.limit=n.limit)
        }
      }else if(method=='rank'){
        ##################################################
        # Rank/TSP Computation of the conditional elements
        if(!is.null(parl)){
          #############################
          #Parallelized VE estimation
          VEs[[j+1]]<-parallel::parApply(clus,indices[[j+1]], 2,estim.VE.rank,
                                         dataX=X,
                                         y=y,
                                         n=n)
          
        }else{
          #############################
          # VE estimation
          VEs[[j+1]]<-apply(indices[[j+1]], 2, estim.VE.rank,
                            dataX=X,
                            y=y,
                            n=n)
        }
      }

    }
  }else{
    #Forces the Shapley effects to sum up to 1
    for(j in 1:(d-1)){
      if(method=='knn'){
        ##################################################
        # KNN Computation of the conditional elements
        if(!is.null(parl)){
          #############################
          #Parallelized VE estimation
          VEs[[j+1]]<-parallel::parApply(clus,indices[[j+1]], 2,estim.VE.knn,
                                         dataX=X,
                                         y=y,
                                         n=n,
                                         n.knn=n.knn,
                                         n.limit=n.limit)
          
        }else{
          #############################
          # VE estimation
          VEs[[j+1]]<-apply(indices[[j+1]], 2, estim.VE.knn,
                            dataX=X,
                            y=y,
                            n=n,
                            n.knn=n.knn,
                            n.limit=n.limit)
        }
    
      }else if (method=="rank"){
        ##################################################
        # Rank/TSP Computation of the conditional elements
        if(!is.null(parl)){
          #############################
          #Parallelized VE estimation
          VEs[[j+1]]<-parallel::parApply(clus,indices[[j+1]], 2,estim.VE.rank,
                                         dataX=X,
                                         y=y,
                                         n=n)
          
        }else{
          #############################
          # VE estimation
          VEs[[j+1]]<-apply(indices[[j+1]], 2, estim.VE.rank,
                            dataX=X,
                            y=y,
                            n=n)
        }
      }
    }  
  VEs[[d+1]]=1
  }
  if(get.sob){
    #Only computing the Sobol' indices
    return(unlist(VEs)[-1])
  }
  ###############################
  #Affecting the weights
  if(!is.null(parl)){
    ################################
    #Parallelized weight affectation
    if(is.null(n.perm)){
      #One core per input
      Shaps <- foreach::foreach(i=1:d, .combine=cbind)%dopar%{
        res=0
        #For every order, the increments are weighted and summed
        #for the input i
        for(ord in 1:d){
          if(ord==1){
            idx<-which(indices[[ord+1]]==i)
            res=res+comb_weights[ord]*VEs[[ord+1]][as.vector(idx)]
          }else{
            #Columns of subsets of order ord+1 containing var_j
            idx_j<-which(indices[[ord+1]]==i, arr.ind=T)[,2]
            #Columns of subsets of order ord without j
            idx_woj<-which(apply(indices[[ord]]!=i, 2, all))
            #Total increment value for order ord
            tot_incr<-sum(VEs[[ord+1]][as.vector(idx_j)]) - sum(VEs[[ord]][as.vector(idx_woj)])
            res=res+comb_weights[ord]*tot_incr
          }
        }
        res
      }
    }else{
      #Parallelized randperm
      Shaps <- foreach::foreach(i=1:d, .combine=cbind)%dopar%{
        res=0
        for(pm in 1:n.perm){
          perm=perms[,pm]
          ord=which(perm==i)
          if(ord!=1){
            idx_j=which(colSums(matrix(apply(indices[[ord+1]], 2 , function(x) perm[1:ord] %in% x), ncol=ncol(indices[[ord+1]])))==ord)
            idx_woj=which(colSums(matrix(apply(indices[[ord]], 2 , function(x) perm[1:(ord-1)] %in% x), ncol=ncol(indices[[ord]])))==ord-1)
            res=res+VEs[[ord+1]][idx_j]-VEs[[ord]][idx_woj]
          }else{
            idx_j=which(indices[[ord+1]]==i)
            res=res+VEs[[ord+1]][idx_j]
          }
        }
        res
      }
    }

  }else{
    ################################
    #Unparallelized weight affectation
    Shaps<-rep(0,d)
    if(is.null(n.perm)){
      for(var_j in 1:d){
        #For every order, the increments are weighted and summed
        #for the input var_j
        for(ord in 1:d){
          if(ord==1){
            idx<-which(indices[[ord+1]]==var_j)
            Shaps[var_j]=Shaps[var_j]+comb_weights[ord]*VEs[[ord+1]][as.vector(idx)]
          }else{
            #Columns of subsets of order ord+1 containing var_j
            idx_j<-which(indices[[ord+1]]==var_j, arr.ind=T)[,2]
            #Columns of subsets of order ord without j
            idx_woj<-which(apply(indices[[ord]]!=var_j, 2, all))
            #Total increment value for order ord
            tot_incr<-sum(VEs[[ord+1]][as.vector(idx_j)]) - sum(VEs[[ord]][as.vector(idx_woj)])
            Shaps[var_j]=Shaps[var_j]+comb_weights[ord]*tot_incr
          }
        }
      }
      
    }else{
      #Unparallelized randperm
      for(var_j in 1:d){
        for(pm in 1:n.perm){
          perm=perms[,pm]
          ord=which(perm==var_j)
          if(ord!=1){
            idx_j=which(colSums(matrix(apply(indices[[ord+1]], 2 , function(x) perm[1:ord] %in% x), ncol=ncol(indices[[ord+1]])))==ord)
            idx_woj=which(colSums(matrix(apply(indices[[ord]], 2 , function(x) perm[1:(ord-1)] %in% x), ncol=ncol(indices[[ord]])))==ord-1)
            Shaps[var_j]=Shaps[var_j]+VEs[[ord+1]][idx_j]-VEs[[ord]][idx_woj]
          }else{
            idx_j=which(indices[[ord+1]]==var_j)
            Shaps[var_j]=Shaps[var_j]+VEs[[ord+1]][idx_j]
          }
        }
      }
    }

  }
  if(is.null(n.perm)){
    Shaps=Shaps/d
  }else{
    Shaps=Shaps/n.perm
  }

  
  if(boot){
    #If Bootstrap computation, return the Shapley Effects
    return(Shaps)
  }else{
    #For non-bootstrap computation, return conditional elements and Shapley Effects
    return(list("VE"=VEs,
                "Shap"=Shaps))
  }
}

shapleysobol_knn <- function(model=NULL, X, method="knn", n.knn = 2, n.limit = 2000, U=NULL, n.perm=NULL, noise=F, rescale=F, nboot = NULL, boot.level = 0.8, conf=0.95, parl=NULL, ...){
  ###########################################
  #Pre-processing
  d<-ncol(X)
  n<-nrow(X)
  
  ###########################################
  #Computing y
  if(!is.null(model)){
    y=model(X,...)
  }else{
    y=NULL
  }
  ##########################################
  #Arguments check
  
  #X
#  if(!(class(X)[1]=="matrix" | class(X)[1]=="data.frame")){
  if(!(inherits(X, "matrix") | inherits(X, "data.frame"))){
    stop("X must be either a matrix of a data frame.")
  }
  #Setting-up colnames
  if(is.null(colnames(X))){
    colnames(X)<-paste("X", 1:d, sep="")
  }
  #Method
  if(!is.element(method, c("rank", "knn"))){
    stop("method must either be rank or knn.")
  }
  
  if(!is.null(n.perm)){
    if(!is.numeric(n.perm)){
      stop("n.perm must be a positive integer.")
    }else if(n.perm>=factorial(d)){
      warning("The number of random permutations is greater than the total number of permutation. Defaulted to n.perm=NULL.")
      n.perm=NULL
    }else if(n.perm<=0){
      stop("n.perm must be a positive integer.")
    }
    else if(n.perm%%1!=0){
      warning("n.perm must be a positive integer. Defaulted to its floor value.")
      n.perm=floor(n.perm)
    }
  }
  #Extracting Sobols
  if(!is.null(U)){
#    if(!class(U)[1]%in%c("numeric","matrix","integer", "list")){
    if(!inherits(U, "numeric") && !inherits(U, "matrix") && !inherits(U, "integer") && typeof(U) != "list"){
      stop("U must be either 0 (Total Sobol' indices), 1 (First-Order indices), a list of subsets, or a matrix.")
#    }else if(class(U)[1]%in%c("integer","numeric")){
    }else if(inherits(U, "integer") || inherits(U,"numeric")){
      if(!U%in%c(0,1)){
        stop("U must be either 0 (Total Sobol' indices), 1 (First-Order indices), a list of subsets, or a matrix.")
      }
    } 
    if(inherits(U, "matrix")){
      if(ncol(U)!=d){
        stop("U must have exactly ncol(X) columns.")
      }
    }
  }
  #Bootstrap Checks
  if(!is.null(nboot)){
    if(!is.numeric(nboot)){
      stop("The 'nboot' argument must be a positive integer.")
    }else if (nboot%%1!=0 | nboot<0 | length(nboot)!=1){
      stop("The 'nboot' argument must be a positive integer.")
    }else if(nboot==0){
      boot=F
    }else{
      boot=T
    }
  }else{
    boot=F
  }

  if(!is.numeric(boot.level)){
    stop("The 'boot.level' argument must be a numeric argument strictly between 0 and 1.")
  }else if(boot.level<=0 | boot.level >=1){
    stop("The 'boot.level' argument must be a numeric argument strictly between 0 and 1.")
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
  
  if(!is.null(y)){
    out<-compute.shapleysobol_knn(X,
                                y,
                                method,
                                n.knn ,
                                n.limit,
                                U,
                                n.perm,
                                noise,
                                rescale,
                                nboot,
                                boot.level,
                                conf,
                                parl,
                                boot)
    out$call=match.call()
  }else{
    out<-list("call"=match.call(),
              "Shap"=NULL,
              "VE"=NULL,
              "indices"=NULL,
              "method"=method,
              "n.perm"=n.perm,
              "w"=NULL,
              "conf_int"=NULL,
              "X"=X,
              "y"=y,
              "n.knn"=n.knn,
              "U"=U,
              "rescale"=rescale,
              "n.limit="=n.limit,
              "boot.level"=boot.level,
              "noise"=noise,
              "boot"=boot,
              "nboot"=nboot,
              "parl"=parl,
              "conf"=conf
    )
    class(out)<-"shapleysobol_knn"
  }
  return(out)
}


compute.shapleysobol_knn<-function(X, y, method, n.knn, n.limit, U, n.perm, noise, rescale, nboot, boot.level, conf, parl, boot){ 
  ###########################################
  #Pre-processing
  d<-ncol(X)
  n<-nrow(X)
  
  ########################################
  #Data Processing
  
  #Identifying categorical inputs
  is.cat<-sapply(1:d, function(x)"factor"%in%class(X[,x]))
  #Rescaling the inputs if specified
  if(rescale & ncol(X[,!is.cat])>1){
    X[,!is.cat]=rescale_inputs(as.matrix(X[,!is.cat]))
  }
  
  #########################################
  #Computing Shapleys
  if(is.null(U)){

  #List of subsets by order ([[2]]: subsets of order 1...etc)
  indices<-rep(list(0), d+1)
  #Weights vector, each element is the weight for |A| ([1]:|A|=0, [2]:|A|=1, ...)
  comb_weights<-rep(0,d)
  
  #Computing subsets and weights
  if(is.null(n.perm)){
    perms=NULL
    ##########################
    #All permutations are used
    for(j in 1:d){
      indices[[j+1]]<-t(gtools::combinations(n=d, r=j))
      comb_weights[j]<-1/choose(n=(d-1), j-1)
    }
  }else{
    comb_weights<-rep(0,d)
    perms=matrix(NA, nrow=d, ncol=n.perm)
    #####################################
    # n.perm random permutations are used
    for(j in 1:d){
      indices[[j+1]]<-matrix(NA, nrow=j, ncol=choose(n=d, j))
      comb_weights[j]<-1
    }
    for(pm in 1:n.perm){
      #We sample a permutation of {1,...,d}
      perm<-sample.int(d,d)
      perms[,pm]<-perm
      #We fill indices accordingly
      for(j in 1:d){
        # subset of order j from the selected random permutation perm
        subs<-perm[1:j]
        if(any(is.na(indices[[j+1]][,1]))){
          #If it's the first subset of order j in indices
          indices[[j+1]][,1]<-subs
        }else{
          chk.sub=colSums(matrix(apply(indices[[j+1]], 2, function(x) subs %in% x), ncol=ncol(indices[[j+1]])))==j
          if(!any(chk.sub)){
              #If the subset is not already in indices, we add it
              id.sub<-max(which(colSums(matrix(is.na(indices[[j+1]]), ncol=ncol(indices[[j+1]])))==0))+1
              indices[[j+1]][,id.sub]<-subs
          }
        }
      }
    }
    #Removing NAs
    for(j in 1:d){
      indices[[j+1]]<-matrix(indices[[j+1]][,stats::complete.cases(t(indices[[j+1]]))], nrow=j)
    }
  }


  
  ########################################
  #Shap KNN Computation
  
  #Preparing parallel cluster
  if(!is.null(parl)){
    cl=parallel::makeCluster(parl)
    doParallel::registerDoParallel(cl)
  }else{
    cl=NULL
  }

  res_Shap<-calc.shapknn(X=X,
                        y=y,
                        method=method,
                        n.perm=n.perm,
                        perms=perms,
                        d=d,
                        indices=indices,
                        comb_weights=comb_weights,
                        parl=parl,
                        clus=cl,
                        boot=F,
                        n.limit=n.limit,
                        n.knn=n.knn,
                        noise=noise)
  
  res_Shap$Shap<-matrix(res_Shap$Shap, ncol=1)
  colnames(res_Shap$Shap)<-"original"
  rownames(res_Shap$Shap)<-colnames(X)
  
  out<-list("call"=match.call(),
            "Shap"=res_Shap$Shap,
            "VE"=res_Shap$VE,
            "indices"=indices,
            "method"=method,
            "n.perm"=n.perm,
            "w"=comb_weights,
            "conf_int"=NULL,
            "X"=X,
            "y"=y,
            "n.knn"=n.knn,
            "rescale"=rescale,
            "n.limit="=n.limit,
            "boot.level"=boot.level,
            "noise"=noise,
            "boot"=boot,
            "nboot"=nboot,
            "parl"=parl,
            "conf"=conf
  )
                           
  if(boot){
    #TODO Mieux travailler le warning.
    warning("Bootstrap confidence intervals may not be theoretically sound.")
    boot.size<-round(n*boot.level)
    sample.boot <- function(data,mle){
      out <- data[sample.int(n=nrow(data), size=mle),,drop=FALSE]
      return(out)
    }
    
    if(!is.null(parl)){
      res.boot <- boot::boot(data=cbind(y,X), statistic=calc.shapknn, 
                             R = nboot, 
                             sim = "parametric", 
                             ran.gen = sample.boot, 
                             mle = boot.size,
                             cl=cl,
                             n.perm=n.perm,
                             perms=perms,
                             ncpus=parl,
                             parallel="snow",
                             y=y,
                             d=d,
                             method=method,
                             indices=indices,
                             comb_weights=comb_weights,
                             parl=NULL,
                             clus=NULL,
                             boot=T,
                             n.limit=n.limit,
                             n.knn=n.knn,
                             noise=noise)
    }else{
      res.boot <- boot::boot(data=cbind(y,X), statistic=calc.shapknn, 
                             R = nboot, 
                             sim = "parametric", 
                             ran.gen = sample.boot, 
                             mle = boot.size,
                             y=y,
                             d=d,
                             method=method,
                             n.perm=n.perm,
                             perms=perms,
                             indices=indices,
                             comb_weights=comb_weights,
                             parl=NULL,
                             clus=NULL,
                             boot=T,
                             n.limit=n.limit,
                             n.knn=n.knn,
                             noise=noise)
      
    }
    
    out$conf_int<-bootstats(res.boot, conf, "basic")
    rownames(out$conf_int)<-colnames(X)
  }
  
  #Stopping the cluster after parallelized computation
  if(!is.null(parl)){parallel::stopCluster(cl)}
  
  #Customize the class of the output for custom methods
  class(out)<-"shapleysobol_knn"
  return(out)
  }else{
    ###########################################
    #Extracting Sobols
    
    #Formatting indices
#    if(class(U)[1]%in%c("integer", "numeric")){
    if(inherits(U, "integer") || inherits(U, "numeric")){
      if(U==0){
        #Only Total Sobol indices
        indices<-rep(list(0),2)
        mat.indices<-matrix(NA, nrow=d-1, ncol=d)
        for(i in 1:d){
          mat.indices[,i]<-setdiff(1:d, i)
        }
        indices[[2]]<-mat.indices
        cname<-"Total Sobol"
      }else if(U==1){
        #Only First-Order Sobol indices
        indices<-rep(list(0),2)
        indices[[2]]<-t(matrix(1:d))
        cname<-"Sobol"
      }
      #Extracting variables names
      rnames<-colnames(X)
    }else{
      #Computing Sobols of subsets
      if(inherits(U, "matrix")){
        #If U is a matrix, transform it to a list of subsets
        U<-apply(U,1,function(x) which(as.logical(x)))
        if(typeof(U)!="list"){
          U<-as.list(U)
        }
      }
      #When U is a list of subsets
      nb.subsets<-length(U)
      #Finding the number of different orders
      ords<-rep(0, nb.subsets)
      for(i in 1:nb.subsets){
        ords[i]<-length(U[[i]])
      }
      un.orders<-unique(ords)
      nb.orders<-length(un.orders)
      indices<-rep(list(0), nb.orders+1)
      #Putting the subsets in indices
      for(i in 1:nb.orders){
        curr.ord<-un.orders[i]
        curr.n.ord<-sum(ords==curr.ord)
        indices[[i+1]]<-matrix(NA, nrow=curr.ord, ncol=curr.n.ord)
        for(j in 1:length(which(ords==curr.ord))){
          sub<-which(ords==curr.ord)[j]
          indices[[i+1]][,j]<-U[[sub]]
        }
      }
      rnames<-unlist(sapply(indices[-1], function(x){
        res<-apply(x, 2, function(y) paste(colnames(X)[y], collapse="-"))
        return(res)
      }))
      cname<-"Closed Sobol"
    }  
    #From indices, extract Sobols
    if(!is.null(parl)){
      #If paralel, setting up the cluster
      cl=parallel::makeCluster(parl)
      doParallel::registerDoParallel(cl)
    }else{
      cl=NULL
    }
    S<-calc.shapknn(X=X, 
                 y=y, 
                 method=method, 
                 n.perm=n.perm, 
                 perms=NULL, 
                 d=d, 
                 indices=indices, 
                 comb_weights=NULL, 
                 parl=parl, 
                 clus=cl, 
                 boot=F, 
                 n.limit=n.limit, 
                 n.knn=n.knn, 
                 noise=T, 
                 i=1:nrow(X), 
                 get.sob=T)
    if(inherits(U, "integer") || inherits(U, "numeric")){
      if(U==0){
        S<-1-S
      }
    }
    S<-matrix(S, ncol=1)
    rownames(S)<-rnames
    colnames(S)<-cname
    out<-list("call"=match.call(),
              "Sobol"=S,
              "indices"=indices[-1],
              "method"=method,
              "conf_int"=NULL,
              "X"=X,
              "y"=y,
              "U"=U,
              "n.knn"=n.knn,
              "rescale"=rescale,
              "n.limit="=n.limit,
              "boot.level"=boot.level,
              "noise"=noise,
              "boot"=boot,
              "nboot"=nboot,
              "parl"=parl,
              "conf"=conf
    )
    class(out)<-"sobol_knn"
    
    if(boot){
      #TODO Mieux travailler le warning.
      warning("Bootstrap confidence intervals may not be theoretically sound.")
      boot.size<-round(n*boot.level)
      sample.boot <- function(data,mle){
        out <- data[sample.int(n=nrow(data), size=mle),,drop=FALSE]
        return(out)
      }
      
      if(!is.null(parl)){
        res.boot <- boot::boot(data=cbind(y,X), statistic=calc.shapknn, 
                               R = nboot, 
                               sim = "parametric", 
                               ran.gen = sample.boot, 
                               mle = boot.size,
                               cl=cl,
                               n.perm=NULL,
                               perms=NULL,
                               ncpus=parl,
                               parallel="snow",
                               y=y,
                               d=d,
                               method=method,
                               indices=indices,
                               comb_weights=NULL,
                               parl=NULL,
                               clus=NULL,
                               boot=T,
                               n.limit=n.limit,
                               n.knn=n.knn,
                               noise=T,
                               get.sob=T)
      }else{
        res.boot <- boot::boot(data=cbind(y,X), statistic=calc.shapknn, 
                               R = nboot, 
                               sim = "parametric", 
                               ran.gen = sample.boot, 
                               mle = boot.size,
                               y=y,
                               d=d,
                               method=method,
                               n.perm=NULL,
                               perms=NULL,
                               indices=indices,
                               comb_weights=NULL,
                               parl=NULL,
                               clus=NULL,
                               boot=T,
                               n.limit=n.limit,
                               n.knn=n.knn,
                               noise=T,
                               get.sob=T)
        
      }
      
      #Bootstraped results for Total indices
      out$conf_int<-bootstats(res.boot, conf, "basic")
      if(inherits(U, "integer") || inherits(U, "numeric")){
        if(U==0){
          out$conf_int[,c(1,4,5)]<-1-out$conf_int[,c(1,4,5)]
          #Switching min/max CI names
          colnames(out$conf_int)<-colnames(out$conf_int)[c(1,2,3,5,4)]
          out$conf_int<-out$conf_int[,c(1,2,3,5,4)]
        }
      }
      rownames(out$conf_int)<-rnames
    }
    #If parallel, stopping the cluster
    if(!is.null(parl)){parallel::stopCluster(cl)}
    return(out)
  }
}


#######################################
# Custom methods

print.shapleysobol_knn<-function(x, ...){
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if(x$method=="knn"){
    cat("\nShapley effects estimation by nearest-neighbor procedure\n")
  }else if(x$method=="rank"){
    cat("\nShapley effects estimation by ranking procedure\n")
  }


  if(x$boot){
    print(x$conf_int)
  }else{
    print(x$Shap)
  }
}


plot.shapleysobol_knn<-function(x, ylim=c(0,1), ...){
  if(x$method=="knn"){
    plot_title="Shapley effects estimation by nearest-neighbor procedure"
  }else if(x$method=="rank"){
    plot_title="Shapley effects estimation by ranking procedure"
  }
  
  plot(x$Shap,
       ylim=ylim,
       axes=F,
       xlab="Inputs",
       ylab="Shapley effects",
       main=plot_title,
       ...)
  axis(2)
  axis(1, at=seq_along(x$Shap), labels=colnames(x$X))
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
           legend=c("Shapley effects",
                    paste(x$conf*100, "% Confidence interval", sep="")),
           bg="white")
  }else{
    legend("topright",
           pch=1,
           legend=c("Shapley effects"),
           bg="white")
  }
}

tell.shapleysobol_knn<-function(x,y,...){
  
  id.obj<-deparse(substitute(x))
  if (! is.null(y)) {
    x$y <- y
  }else if (is.null(x$y)) {
    stop("y not found.")
  }
  
  res<-compute.shapleysobol_knn(X=x$X,
                              y=x$y,
                              method=x$method,
                              n.knn=x$n.knn ,
                              n.limit=x$n.limit,
                              U=x$U,
                              n.perm=x$n.perm,
                              noise=x$noise,
                              rescale=x$rescale,
                              nboot=x$nboot,
                              boot.level=x$boot.level,
                              conf=x$conf,
                              parl=x$parl,
                              boot=x$boot)
  assign(id.obj, res, parent.frame())
  
  return(res)
}

extract.shapleysobol_knn<- function(x, ...){
  if (is.null(x$VE) | is.null(x$Shap)){
    stop("The Shapley effects must be computed (U=NULL) prior to extracting first-order and total Sobol' indices.")
  }
  d<-ncol(x$X)
  Sob<-as.matrix(x$VE[[2]], ncol=1)
  Sob.tot<-1-as.matrix(rev(x$VE[[d]]), ncol=1)
  colnames(Sob)<-"Sobol"
  rownames(Sob)<-rownames(x$Shap)
  colnames(Sob.tot)<-"Total Sobol"
  rownames(Sob.tot)<-rownames(x$Shap)
  return(list("Sobol"=Sob,
              "Sobol.tot"=Sob.tot))
}

#Plotting the Sobol' indices
plot.sobol_knn<-function(x, ylim=c(0,1), ...){
  if(inherits(x$U, "integer") || inherits(x$U, "numeric")){
    if(x$U==0){
      if(x$method=="knn"){
        plot_title=("\nTotal Sobol' indices estimation by nearest-neighbor procedure\n")
      }else if(x$method=="rank"){
        plot_title=("\nTotal Sobol' indices estimation by ranking procedure\n")
      }
      
    }else if(x$U==1){
      if(x$method=="knn"){
        plot_title=("\nFirst order Sobol' indices estimation by nearest-neighbor procedure\n")
      }else if(x$method=="rank"){
        plot_title=("\nFirst order Sobol' indices estimation by ranking procedure\n")
      }
    }
  }else{
    if(x$method=="knn"){
      plot_title=("\nClosed Sobol' indices estimation by nearest-neighbor procedure\n")
    }else if(x$method=="rank"){
      plot_title=("\nClosed Sobol' indices estimation by ranking procedure\n")
    }
  }
  
  plot(x$Sobol,
       ylim=ylim,
       axes=F,
       xlab="Inputs/Subsets",
       ylab="Sobol' indices",
       main=plot_title,
       ...)
  axis(2)
  axis(1, at=seq_along(x$Sobol), labels=rownames(x$Sobol))
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
           legend=c("Sobol' indices",
                    paste(x$conf*100, "% Confidence interval", sep="")),
           bg="white")
  }else{
    legend("topright",
           pch=1,
           legend=c("Sobol' indices"),
           bg="white")
  }
}

#Printing the Sobol' indices
print.sobol_knn<-function(x, ...){
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if(inherits(x$U, "integer") || inherits(x$U, "numeric")){
    if(x$U==0){
      if(x$method=="knn"){
        cat("\nTotal Sobol' indices estimation by nearest-neighbor procedure\n")
      }else if(x$method=="rank"){
        cat("\nTotal Sobol' indices estimation by ranking procedure\n")
      }
      
    }else if(x$U==1){
      if(x$method=="knn"){
        cat("\nFirst order Sobol' indices estimation by nearest-neighbor procedure\n")
      }else if(x$method=="rank"){
        cat("\nFirst order Sobol' indices estimation by ranking procedure\n")
      }
    }
  }else{
    if(x$method=="knn"){
      cat("\nClosed Sobol' indices estimation by nearest-neighbor procedure\n")
    }else if(x$method=="rank"){
      cat("\nClosed Sobol' indices estimation by ranking procedure\n")
    }
  }
  
  if(x$boot){
    print(x$conf_int)
  }else{
    print(x$Sobol)
  }
}

ggplot.shapleysobol_knn<-function(x, ylim = c(0, 1), ...) {
  if (!is.null(x$Shap)){
    shap <- as.data.frame(x$Shap) 
    colnames(shap) <- "original"
    nodeggplot(list(shap), xname = "Shapley effects", ylim = ylim)
  }else{
    stop("No plot available")
  }
}

ggplot.sobol_knn<-function(x, ylim=c(0,1), ...){
  if(!is.null(x$Sobol)){
    sobols<-as.data.frame(x$Sobol)
    colnames(sobols)<-"original"
    if(typeof(x$U)=="list"){
      x_name<-"Closed Sobol'"
    }else{
      if(x$U==1){
        x_name<-"First Order Sobol'"
      }else if (x$U==0){
        x_name<-"Total Sobol'"
      }
    }
    nodeggplot(list(sobols), xname=x_name, ylim=ylim)
  }else{
    stop("No plot available")
  }

}