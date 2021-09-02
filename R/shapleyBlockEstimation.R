
shapleyBlockEstimationS <- function(Beta, S, kappa=0,  M=20, tol=10^(-6))
{
  #####################################################
  # This function computes the Shapley effects in the linear Gaussian framework 
  # estimating the block-diagonal structure of the covariance matrix of the inputs
  #
  # List of inputs to this function:
  # Beta: a vector containing the coefficients of the linear model (without the value at zero).
  # S: the empirical covariance matrix of the inputs. Has to be positive semi-definite matrix with same size that Beta.
  # kappa: the penalization coefficient that promotes block-diagonal structures.
  # M: maximal size of the estimate of the block-diagonal structure.
  # tol: a relative tolerance to detect zero singular values of Sigma.
  #
  #####################################################
  
  estim=EstimCov(S, kappa, M=M) 
  eta=shapleyGaussianGroups(Beta, S, estim$label,tol=tol)
  estim$Shapley=eta
  
  return(estim)
}


shapleyBlockEstimationX <- function(X,Y,delta=NULL, M=20, tol=10^(-6))
{
  #####################################################
  # This function computes the Shapley effects in the linear Gaussian framework 
  # estimating the block-diagonal structure of the covariance matrix of the inputs
  #
  # List of inputs to this function:
  # X: a matrix containing an i.i.d. sample of the inputs.
  # Y: a vector containing the corresponding i.i.d. sample of the (noisy) output.
  # delta: a positive number that fixes the positive penalization coefficient to 1/(p n^delta)
  # M: maximal size of the estimate of the block-diagonal structure.
  # tol: a relative tolerance to detect zero singular values of Sigma.
  #
  #####################################################
  X=as.matrix(X)
  dimen=dim(X)
  n=dimen[1]
  p=dimen[2]
  if(n<=p)
  {
    stop('error: "n" is not greater than "p"')
  }
  if(length(Y)!=n)
  {
    stop('error: "X" and "Y" do not have good dimensions')
  }
  
  S=var(X)
  
  # estimation of the linear model
  XX=cbind(rep(1,n),X)
  beta=solve(t(XX) %*% XX,t(XX)%*%as.matrix(Y))
  Beta=beta[-1]
  
  # choice of the penalization coefficient using the value of delta
  if(is.null(delta)){
    kappa=0
  }else{
    kappa=1/(p*n^(delta))
  }  
  
  estim=shapleyBlockEstimationS(Beta=Beta,S=S,kappa=kappa,M=M,tol=tol)
  return(estim)
}




shapleyGaussianGroups <-function(Beta,Sigma,groups,tol=10^(-6)){
  
  K=max(groups)
  
  p=length(Beta)
  eta=rep(0,p)
  
  VarY=t(Beta)%*%Sigma%*%Beta
  
  Eta=function(Beta_k,Sigma_k)
  {
    p_k=length(Beta_k)
    if(p_k==1)
    {
      return(Beta_k^2*Sigma_k)
    }
    
    Eu=rep(0,2^p_k)
    Eu[2^p_k]=t(Beta_k)%*%Sigma_k%*%Beta_k
    for(j in 1:(2^p_k-2))
    {
      u=as.numeric(intToBits(j)[1:p_k])
      Eu[j+1]=Beta_k[u==1]%*%(Sigma_k[u==1,u==1]-(Sigma_k[u==1,u==0]%*%
                                                    MASS::ginv(Sigma_k[u==0,u==0],tol=tol)%*%Sigma_k[u==0,u==1]))%*%Beta_k[u==1]
    }
    eta_k=rep(0,p_k)
    for (i in 1:p_k)
    {
      for(l in 0: (2^p_k-1))
      {
        if (floor(l/(2^(i-1)))%%2==0)
        {
          u=as.numeric(intToBits(l)[1:p_k])
          cardu=sum(u)
          eta_k[i]=eta_k[i]+(Eu[l+1+2^(i-1)]- Eu[l+1])/choose(p_k-1,cardu)
        }
      }
    }
    eta_k=eta_k/p_k
    return(eta_k)
  }
  
  
  eta=rep(0,p)
  for (k in 1:K)
  {
    Ind=which(groups==k)
    eta_k=Eta(Beta[Ind],Sigma[Ind,Ind])
    eta[Ind]=eta_k
  }
  eta=as.vector(eta)/as.vector(VarY)
  return(eta)
}




EstimCov=function(S,kappa,M=NULL,cond=FALSE)
{
  p=dim(S)[1]
  V=diag(1/sqrt(diag(S)))
  cor=abs(V%*%S%*%V) #valeur absolue de la matrice de correlation
  cor[lower.tri(cor)]=t(cor)[lower.tri(cor)] # on rend la matrice symmetrique (elle ne l'etait pas a cause des erreurs)
  
  valThres <- cor[upper.tri(cor)] 
  thres <- rev(valThres[order(valThres)])
  
  
  labelsPath <- list()
  
  
  lam=1
  len=length(thres)
  repeat{
    lambdaR <- thres[lam]
    E <- cor
    E[cor>lambdaR] <- 1
    E[cor<lambdaR] <- 0
    E[cor==lambdaR] <- 0
    goutput <- igraph::graph.adjacency(E,mode="undirected",weighted=NULL)
    label= igraph::clusters(goutput)$membership
    m=max(table(label))
    if (m>M || lam==len)
      break
    labelsPath[[lam]] <- label
    lam=lam+1
  }   
  
  
  
  labels=unique(labelsPath) # we remove the groups that occur several times
  
  G=length(labels)
  if(kappa==0)
  {
    label=labels[[G]]
  }else{
    Max=unlist(lapply(labels, function(x) max(table(x))))
    
    Phi<-function (groups) # -log likelihood + penalization
    {
      K=max(groups)
      phi=0
      for(k in 1:K)
      {
        Ind=which(groups==k)
        phi=phi+log(det(as.matrix(S[Ind,Ind])))/p+kappa*length(Ind)^2
      }
      return(phi)
    }
    
    PHI=unlist(lapply(labels,Phi))
    
    Estim=list("labels"=labels, "Max"=Max, "Phi"=PHI)
    
    l=which.min(Estim$Phi)
    label=Estim$labels[[l]]
  }
  K=max(label)
  
  S_B=matrix(0,nrow=p,ncol=p)
  for(k in 1:K)
  {
    Ind=which(label==k)
    S_B[Ind,Ind]=S[Ind,Ind]
  }
  out=list("label"=label,"S_B"=S_B)
  
  return(out)
}



