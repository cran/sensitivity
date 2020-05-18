shapleyLinearGaussian <-function(Beta,Sigma,tol=10^(-6)){
  #####################################################
  # This function computes the Shapley effects in the linear Gaussian framework
  #
  # List of inputs to this function:
  # Beta: a vector containing the coefficients of the linear model (without the value at zero).
  # Sigma: covariance matrix of the inputs. Has to be positive semi-definite matrix with same size that Beta. 
  # tol: a relative tolerance to detect zero singular values of Sigma.
  #
  #####################################################
  
  
  # require
  # library(MASS)
  # library(igraph)
  
  
  
  p=length(Beta)
  mn=dim(Sigma)
  if(mn[1]!=p | mn[2]!=p)
  {
    print("'Sigma' should be a symmetric matrix with same size than 'Beta'.")
  }

  
  VarY=t(Beta)%*%Sigma%*%Beta
  
  G=igraph::clusters(igraph::graph.adjacency(Sigma,weighted = TRUE))
  
  if(max(G$csize)>=25)
  {
    rep=readline("The largest block has a size larger than 25. Are you sure to continue? (y/n) ")
    if(rep!="y")
    {
      return(0)
    }
  }
  
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
  K=G$no
  for (k in 1:K)
  {
    eta_k=Eta(Beta[G$membership==k],Sigma[G$membership==k,G$membership==k])
    eta[G$membership==k]=eta_k
  }
  eta=as.vector(eta)/as.vector(VarY)
  return(eta)
}