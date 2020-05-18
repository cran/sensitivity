
shapleySubsetMc <- function(X, Y, Ntot=NULL, Ni=3, cat=NULL, weight=NULL, discrete=NULL, noise=FALSE) {
  #####################################################
  # This function estimates the Shapley effects from data
  #
  # List of inputs to this function:
  # X: sample of inputs
  # Y: sample of output
  # Ntot: around sum of the costs Nu of the estimates of the Vu
  #       if Ntot=NULL, the maximum cost is choosen (advised for little datas)
  # Ni: accuracy of the first loop Monte-Carlo
  # cat: indicates the categorical inputs (factor variables, not real variables).
  #      the real variables do not need to be continuous
  # weight: if weight=NULL, all the categorical inputs will have the same weight.
  #         Otherwise, weight is a vector with the same length of cat and with 
  #         the weight of each input variable.
  # discrete: indicate the real variables that give a positive probability to some points.
  #         the categorical variables do not have to be in 'discrete'
  # noise: indicates if "Y=f(X)" or if "Y=f(X)+ noise" and so if E(V(Y|X)) needs to be estimated
  #
  #####################################################


  # First Step: we estimate the Vu

  #initialization

  dimen=dim(X)
  N=dimen[1]
  p=dimen[2]
  if(length(Y)!=N)
  {
    stop("the inputs and the output do not have the same number of individuals")
  }
  if(length(weight)!=0 & length(weight)!=length(cat))
  {
    stop("'weight' must have the same length as 'cat'")
  }


  if(mode(X)!="numeric")
  {
    if(mode(X)=="list")
    {
      X=as.matrix(sapply(X, as.numeric))
    }else{
      stop("'X' is not a matrix or a dataframe")
    }
  }

  if(is.null(colnames(X)))
  {
    names=paste(rep("X",p),1:p,sep="")
  }else{
    names=colnames(X)
  }

  # we add the categorical variables in the discrete variables
  discr=union(cat,discrete)

  EY <- mean(Y)
  VarY=var(Y)

  Varu=rep(0,2^p)
  


  U=t(sapply(1:(2^p-2),function(x){as.numeric(intToBits(x)[1:p])}))


  # distance in a product space

  if (length(weight)==0)
  {
    weight=rep(1,length(cat))
  }

  norm_vec <- function(x1,x2,realvar,catvar)
  {
    dist_prod=(sum((x1[realvar]-x2[realvar])^2))+sum((weight[is.element(cat,catvar)])*(x1[catvar]!=x2[catvar])) # we can add sqrt behind the first parenthesis

    return(dist_prod)
  }


  # function which estimates each Vu
  VU=function(u)
  {
    cardu=sum(u)
    if(length(Ntot)==0)
    {
      Nu=N
    }else{
      Nu=round(Ntot/choose(p,cardu)/(p-1)) #accuracy of Vu
    }
    
    if(Nu==0)
    {
      warning("Ntot is too small and some conditional elements have been estimated to 0")
      return(0)
    }else{
      
      # uniform indices of Xu
      Num=sample(1:N,Nu,replace=FALSE)
      
      if(Nu==1)
      {
        X_unif=matrix(X[Num,],nrow=1)
      }else{
        X_unif=X[Num,]
      }
      
      Su=which(u==1) # set u
      Smu=which(u==0) # set -u
      
      cat_Smu=intersect(Smu,cat)
      real_Smu=setdiff(Smu,cat)
      
      V=NULL
      
      # TF is FALSE if all the variables -u are discrete
      TF=(length(setdiff(Smu,discr))>0)
      
      if (TF)
      {
        for (n in 1:Nu)
        {
          xx=X_unif[n,]
          
          dist=apply(X,1,function(x) norm_vec(x,xx,real_Smu,cat_Smu))
          num_cond=order(dist,decreasing = F)[1:Ni]
          
          Vn=var(Y[num_cond])
          V=c(V,Vn)
        }
      }else{
        for (n in 1:Nu)
        {
          xx=X_unif[n,]
          
          # we mix the individuals to avoid taking the same nearest neigbours
          # when the variables -u are all discrete
          samp=sample(1:N)
          XX=X[samp,]
          YY=Y[samp]
          
          
          
          dist=apply(XX,1,function(x) norm_vec(x,xx,real_Smu,cat_Smu))
          num_cond=order(dist,decreasing = F)[1:Ni]
          
          Vn=var(YY[num_cond])
          V=c(V,Vn)
        }
      }
      return(mean(V))
    }
  }
  Varu[2:(2^p-1)]=apply(U,1,VU)
  
  Varu[2^p]=VarY
  
  if(noise)
  {
    Varu[1]=VU(rep(0,p))
    VarY=VarY-Varu[1]
  }
  


  # we compute the real total cost (for information)
  if(length(Ntot)==0)
  {
    if(noise)
    {
      cost=(2^p-1)*N
    }else{
      cost=(2^p-2)*N
    }
  }else{
    cost=sum(apply(U,1,function(u) round(Ntot/choose(p,sum(u))/(p-1))))
    if(noise)
    {
      cost=cost+ round(Ntot/choose(p,0)/(p-1))
    }
  }


  # Second step: we deduce the estimates of the Shapley effects

  eta=rep(0,p)

  for (i in 1:p)
  {
    funci=function(k)
    {
      if (floor(k/(2^(i-1)))%%2==0) # if i is not in u
      {
        return((Varu[k+1+2^(i-1)]-Varu[k+1])/choose(p-1,sum(U[k,])))
      }else{
        return(0)
      }
    }

    diffVa=apply(matrix(0:(2^p-1),nrow=1),2,funci)
    eta[i]=sum(diffVa)/p/VarY
  }

  Shap=list("shapley"=eta, "cost"=cost, "names"=names)
  class(Shap) <- "shapleySubsetMc"
  return(Shap)
}

plot.shapleySubsetMc <- function(x, ylim = c(0, 1), ...) {
  if (!is.null(x$shapley)& !is.null(x$names)) {
    p=length(x$shapley)
    pch = 21
    plot(x$shapley, ylim = ylim, pch = pch,xaxt="n",
         xlab = "", ylab = "")

    axis(side = 1, at=1:p,labels =x$names)
  }
}

