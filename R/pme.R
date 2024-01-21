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

#Returns the Total Sobol' indices from the closed Sobol' indices
get.SobolT<-function(VEs){
  d<-length(VEs)-1
  SobT<-rep(list(0), length(VEs))
  for (ord in 1:d){
    SobT[[ord+1]] = rev(VEs[[d+1]] - VEs[[d+1-ord]])
  }
  return(SobT)
}


#Recursive PV estimation
recur.PV<-function(indices, VEs, parl, cl){
  
  if(!is.null(parl)){doParallel::registerDoParallel(cl)}
  d=length(indices)-1
  Ps<-rep(list(0), d)
  Ps[[1]]=VEs[[2]]
  for(ord in 2:d){
    Ws<-VEs[[ord+1]]
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
  
  return(Ps)
}


########################
#Estimation of all the PVs

calc.pv<-function(X, y, method, d, indices, parl, clus, boot, n.limit, n.knn, noise, tol, marg, i=1:nrow(X)){
                           
  ###############################
  #Pre-Processing
  #Row selection (for Bootstrap estimate, default doesn't change anything)
  if(boot){
    y<-X[,1]
    X<-X[,-1]
  }
  n<-nrow(X)
  
  #List of VEs, arranged by orders (same structure as the indices object)
  VEs<-rep(list(0), d+1)
  
  ####################################
  #Conditional elements estimation
  #For every order, the corresponding VE estimate of each
  #subset of indices are stored in VEs
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
  #Allowing the effects not to sum up to one, but to V(E[Y|X])
  if(noise){
    #If noise == T, then we estimate V(E[Y|X])
    if(method=="knn"){
      VEX<-estim.VE.knn(dataX=X,
                        y=y,
                        subset=indices[[d+1]],
                        n=n,
                        n.knn=n.knn,
                        n.limit=n.limit)
      }else if (method=="rank"){
        VEX<-estim.VE.rank(dataX=X,
                          y=y,
                          subset=indices[[d+1]],
                          n=n)
      }
      VEs[[d+1]]=VEX
    }else{
    #If noise == F, V(E[Y|X]) is equal to 1
    VEs[[d+1]]=1
  }
  
  #If PME computation, then change the closed Sobol' indices to the Total Sobol' indices (their marginal)
  if(marg){
    VEs<-get.SobolT(VEs)
  }
  #################################################
  #PV Computation
  
  #Identifying inputs with 0 Total Sobol (zero players)
  if(!is.null(tol)){
    idx.z<-indices[[2]][,which(abs(VEs[[2]])<=tol)]
  }else{
    idx.z<-indices[[2]][,which(VEs[[2]]==0)]
  }
  if(any(idx.z)){
    #Checking if there are other zero (or under tol) coalitions in each order
    if(!is.null(tol)){
      Z.coal.order<-which(sapply(VEs, function(x) any(abs(x)<=tol)))
      #Which coalitions are null coalitions, may be useful someday
      #Z.coal.idx<-sapply(a.coal.order, function(x) return(indices[[x]][,which(VEs[[x]]<=tol),drop=FALSE]), simplify=F)
    }else{
      Z.coal.order<-which(sapply(VEs, function(x) any(x==0)))
      #Which coalitions are null coalitions, may be useful someday
      #Z.coal.idx<-sapply(a.coal.order, function(x) return(indices[[x]][,which(VEs[[x]]<=tol),drop=FALSE]), simplify=F)
    }
    #Getting the cardinality of the biggest zero (or under tol) coalition
    Z.cardMax<-max(Z.coal.order)-1 #It should always be over or equal to 1 since there are zero players
    #Biggest zero coalitions of max cardinality
    if(!is.null(tol)){
      Z.coal.cardMax<-indices[[Z.cardMax+1]][,which(abs(VEs[[Z.cardMax+1]])<=tol), drop=F]
    }else{
      Z.coal.cardMax<-indices[[Z.cardMax+1]][,which(VEs[[Z.cardMax+1]]==0), drop=F]
    }
    #Boolean of size card(idx.z) indicating if they are in all the zero coalitions of max cardinality
    z.zeroPV<-sapply(idx.z, function(x) sum(colSums(Z.coal.cardMax==x))==1)
    #First case: Z.cardMax=d, no one gets anything (should not happen too often)
    if(Z.cardMax==d){
      PV=rep(0,d)
      PV=matrix(PV, ncol=1)
    }else if(Z.cardMax==d-1){
      #Second case: Z.cardMax=d-1, no need to loop
      PV=rep(0,d)
      idx.pv=setdiff(1:d,idx.z[z.zeroPV])
      PV[idx.pv] = VEs[[d+1]]/nrow(Z.coal.cardMax)
      PV=matrix(PV, ncol=1)
    }else{
      #Third case: Recursive computing
      PS.i<-rep(0,d)
      PS.N<-rep(0,d)
      #Setting up VEs and indices by removing every zero coalition
      for(idx.Zcoal in 1:ncol(Z.coal.cardMax)){
        Z.coal=Z.coal.cardMax[,idx.Zcoal]
        indices_<-rep(list(0), d-Z.cardMax+1)
        VEs_<-rep(list(0), d-Z.cardMax+1)
        for(i in 1:(d-Z.cardMax)){
          checkmat<-matrix(is.element(indices[[i+1]], Z.coal), i, ncol(indices[[i+1]]))
          idx_ind_null <- as.vector(which(colSums(checkmat) == 0))
          ind_tmp<-indices[[i+1]][, idx_ind_null, drop=F]
          indices_[[i+1]]<-rbind(ind_tmp, matrix(rep(Z.coal, ncol(ind_tmp)), ncol=ncol(ind_tmp), nrow=length(Z.coal)))
        }
        for(i in 1:(length(indices_)-1)){
          idx.get<-apply(indices_[[i+1]],
                         2,
                         function(x) which(colSums(matrix(indices[[Z.cardMax+i+1]] %in% x,
                                                          nrow=Z.cardMax+i, 
                                                          ncol=ncol(indices[[Z.cardMax+i+1]])))==i+Z.cardMax)
                         )
          
          VEs_[[i+1]]<-VEs[[Z.cardMax+i+1]][idx.get]
        }
        PS<-recur.PV(indices=indices_,
                     VEs=VEs_, 
                     parl=parl,
                     cl=clus)
        idx.var<-as.vector(indices_[[2]][1,])
        PS.i[idx.var]=PS.i[idx.var]+(1/rev(PS[[length(PS)-1]]))
        PS.N=PS.N+1/PS[[length(PS)]]
      }
      PV=matrix(PS.i/PS.N, ncol=1)
    }

  }else{
    PS<-recur.PV(indices=indices,
                 VEs=VEs,
                 parl=parl, 
                 cl=clus)
    PV<-matrix(rev(PS[[d]]/PS[[d-1]]), ncol=1)
  }
  
  if(!boot){
    return(list("PME"=PV,
                "VE"=VEs))
  }else{
    return(as.vector(PV))
  }

}


pme_knn <- function(model=NULL, X, method="knn", tol=NULL, marg=T, n.knn = 2, n.limit = 2000, noise=F, rescale=F, nboot = NULL, boot.level = 0.8, conf=0.95, parl=NULL, ...){
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
  if(!(class(X)[1]=="matrix" | class(X)[1]=="data.frame")){
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
    out<-compute.pme_knn(X=X,
                         y=y,
                         method=method,
                         tol=tol,
                         marg=marg,
                         n.knn=n.knn,
                         n.limit=n.limit,
                         noise=noise,
                         rescale=rescale,
                         nboot=nboot,
                         boot.level=boot.level,
                         conf=conf,
                         parl=parl,
                         boot=boot)
    
    out$call=match.call()
  }else{
    out<-list("call"=match.call(),
              "PME"=NULL,
              "VE"=NULL,
              "indices"=NULL,
              "method"=method,
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
              "conf"=conf,
              "marg"=marg,
              "tol="=tol
    )
    class(out)<-"pme_knn"
  }
  return(out)
}



compute.pme_knn<-function(X, y, method, tol, marg, n.knn, n.limit, U, noise, rescale, nboot, boot.level, conf, parl, boot){ 
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
  #Computing PME
  #List of subsets by order ([[2]]: subsets of order 1...etc)
  indices<-rep(list(0), d+1)
  #Weights vector, each element is the weight for |A| ([1]:|A|=0, [2]:|A|=1, ...)
    
  #Computing subsets and weights
  ##########################
  #All permutations are used
  for(j in 1:d){
    indices[[j+1]]<-t(gtools::combinations(n=d, r=j))
  }

  ########################################
  #PME KNN Computation
    
  #Preparing parallel cluster
  if(!is.null(parl)){
    cl=parallel::makeCluster(parl)
    doParallel::registerDoParallel(cl)
  }else{
    cl=NULL
  }
    
  res_PME<-calc.pv(X=X,
                      y=y,
                      method=method,
                      d=d,
                      indices=indices,
                      parl=parl,
                      clus=cl,
                      boot=F,
                      n.limit=n.limit,
                      n.knn=n.knn,
                      noise=noise,
                      marg=marg,
                   tol=tol)
    
    colnames(res_PME$PME)<-"original"
    rownames(res_PME$PME)<-colnames(X)
    
    out<-list("call"=match.call(),
              "PME"=res_PME$PME,
              "VE"=res_PME$VE,
              "indices"=indices,
              "method"=method,
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
              "conf"=conf,
              "marg"=marg,
              "tol"=tol
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
        res.boot <- boot::boot(data=cbind(y,X), statistic=calc.pv, 
                               R = nboot, 
                               sim = "parametric", 
                               ran.gen = sample.boot, 
                               mle = boot.size,
                               cl=cl,
                               ncpus=parl,
                               parallel="snow",
                               y=y,
                               d=d,
                               method=method,
                               indices=indices,
                               parl=NULL,
                               clus=NULL,
                               boot=T,
                               n.limit=n.limit,
                               n.knn=n.knn,
                               noise=noise,
                               marg=marg,
                               tol=tol)
      }else{
        res.boot <- boot::boot(data=cbind(y,X), statistic=calc.pv, 
                               R = nboot, 
                               sim = "parametric", 
                               ran.gen = sample.boot, 
                               mle = boot.size,
                               y=y,
                               d=d,
                               method=method,
                               indices=indices,
                               parl=NULL,
                               clus=NULL,
                               boot=T,
                               n.limit=n.limit,
                               n.knn=n.knn,
                               noise=noise,
                               marg=marg,
                               tol=tol)
      }
      
      out$conf_int<-bootstats(res.boot, conf, "basic")
      rownames(out$conf_int)<-colnames(X)
    }
    
    #Stopping the cluster after parallelized computation
    if(!is.null(parl)){parallel::stopCluster(cl)}
    
    #Customize the class of the output for custom methods
    class(out)<-"pme_knn"
    return(out)
}


#######################################
# Custom methods

print.pme_knn<-function(x, ...){
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if(x$method=="knn"){
    cat("\nProportional marginal effects estimation by nearest-neighbor procedure\n")
  }else if(x$method=="rank"){
    cat("\nProportional marginal effects estimation by ranking procedure\n")
  }
  
  if(x$boot){
    print(x$conf_int)
  }else{
    print(x$PME)
  }
}


plot.pme_knn<-function(x, ylim=c(0,1), ...){
  if(x$method=="knn"){
    plot_title="Proportional marginal effects estimation by nearest-neighbor procedure"
  }else if(x$method=="rank"){
    plot_title="Proportional marginal effects estimation by ranking procedure"
  }
  
  plot(x$PME,
       ylim=ylim,
       axes=F,
       xlab="Inputs",
       ylab="Proportional marginal effects",
       main=plot_title,
       ...)
  axis(2)
  axis(1, at=seq_along(x$PME), labels=colnames(x$X))
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
           legend=c("Proportional marginal effects",
                    paste(x$conf*100, "% Confidence interval", sep="")),
           bg="white")
  }else{
    legend("topright",
           pch=1,
           legend=c("Proportional marginal effects"),
           bg="white")
  }
}

ggplot.pme_knn<-function(data, mapping = aes(), ylim = c(0,1), ..., environment = parent.frame()) {
  x <- data
  if (!is.null(x$PME)){
    pme <- as.data.frame(x$PME) 
    colnames(pme) <- "original"
    nodeggplot(list(pme), xname = "Proportional marginal effects", ylim = ylim)
  }else{
    stop("No plot available")
  }
}

tell.pme_knn<-function(x,y,...){
  
  id.obj<-deparse(substitute(x))
  if (! is.null(y)) {
    x$y <- y
  }else if (is.null(x$y)) {
    stop("y not found.")
  }
  
  res<-compute.pme_knn(X=x$X,
                                y=x$y,
                                method=x$method,
                                n.knn=x$n.knn ,
                                n.limit=x$n.limit,
                                noise=x$noise,
                                rescale=x$rescale,
                                nboot=x$nboot,
                                boot.level=x$boot.level,
                                conf=x$conf,
                                parl=x$parl,
                                boot=x$boot,
                      marg=x$marg,
                      tol=x$tol)
  assign(id.obj, res, parent.frame())
  
  return(res)
}

