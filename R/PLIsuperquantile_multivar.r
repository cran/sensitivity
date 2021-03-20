# library(evd)

PLIsuperquantile_multivar = function(order,x,y,inputs,deltasvector,InputDistributions,samedelta=TRUE,percentage=TRUE,nboot=0,conf=0.95,bootsample=TRUE,bias=TRUE){
  
  # This function allows the estimation of Density Modification Based Reliability Sensitivity Indices
  # called PLI (Perturbed-Law based sensitivity Indices) for a superquantile
  #
  # It is adapted for taking into account simultaneous perturbations of 2 inputs
  #
  # It only works for MEAN perturbation of the inputs
  #
  ###################################
  ## Description of input parameters
  ###################################
  #  
  # order is the order of the superquantile to estimate.
  # x is the matrix of simulation points coordinates (one column per variable).
  # y is the vector of model outputs.
  # inputs is the vector of the two inputs' indices
  # deltasvector is a vector containing the values of delta for which the indices will be computed
  # InputDistributions is a list of list. Each list contains, as a list, the name of the distribution to be used and the parameters.
  #   Implemented cases so far: Gaussian, Uniform, Triangle, Left Trucated Gaussian
  # samedelta is a boolean used with the value "MOY" for type. If it is set at TRUE, the mean perturbation will be the same for all the variables. 
  #   If not, the mean perturbation will be new_mean = mean+sigma*delta where mean, sigma are parameters defined in InputDistributions and delta is a value of deltasvector. See subsection 4.3.3 of the reference for an exemple of use. 
  # percentage defines the formula used for the PLI. If percentage=FALSE, the classical formula used in the bibliographic references is used.
  #   If percentage=TRUE, the PLI is given in percentage of variation of the superquantile (even if it is negative). 
  # nboot is the number of bootstrap replicates
  # conf is the required bootstrap confidence interval
  # bootsample defines if the sampling uncertainty is taken into account in computing the boostrap confidence intervals of the PLI
  # bias defines which type of PLI-superquantile is computed:
  #   "TRUE" gives the mean of outputs above the perturbed quantile (alternative formula)
  #   "FALSE" gives the the mean of perturbed outputs above the perturbed quantile (original formula)
  #
  ###################################
  ## Description of output parameters
  ###################################
  #
  # The output is a list of matrices where the PLI and superquantiles are stored
  #
  ########################################
  ## Creation of local variables	
  ########################################
  
  nmbredevariables=dim(x)[2]		# number of variables
  nmbredepoints=dim(x)[1]			# number of points
  nmbrededeltas=length(deltasvector)	# number of perturbations
  
  ## some storage matrices 
  I <- J <- ICIinf <- ICIsup <- JCIinf <- JCIsup <- matrix(0,ncol=nmbrededeltas,nrow=nmbrededeltas) 

  ########################################
  ## Definition of useful functions used further on 
  ########################################
  
  ########################################
  ##	Simpson's method
  ########################################
  simpson_v2 <- function(fun, a, b, n=100) {
    # numerical integral using Simpson's rule
    # assume a < b and n is an even positive integer
    if (a == -Inf & b == Inf) {
      f <- function(t) (fun((1-t)/t) + fun((t-1)/t))/t^2
      s <- simpson_v2(f, 0, 1, n)
    } else if (a == -Inf & b != Inf) {
      f <- function(t) fun(b-(1-t)/t)/t^2
      s <- simpson_v2(f, 0, 1, n)
    } else if (a != -Inf & b == Inf) {
      f <- function(t) fun(a+(1-t)/t)/t^2
      s <- simpson_v2(f, 0, 1, n)
    } else {
      h <- (b-a)/n
      x <- seq(a, b, by=h)
      y <- fun(x)
      y[is.nan(y)]=0
      s <- y[1] + y[n+1] + 2*sum(y[seq(2,n,by=2)]) + 4 *sum(y[seq(3,n-1, by=2)])
      s <- s*h/3
    }
    return(s)
  }
  
  transinverse <- function(a,b,c) {
    # useful to compute bootstrap-based confidence intervals
    if (a > c) ans <- a / b - 1
    else ans <- 1 - b / a
    return(ans)
  }
  
  ########################################
  ## Principal loop of the function 
  ########################################
  
  quantilehat <- quantile(y,order) # quantile estimate
  sqhat <- mean( y[y >= quantilehat] ) # superquantile estimate
  
  ys = sort(y,index.return=T) # ordered output
  xs = x[ys$ix,] # inputs ordered by increasing output
  
  lqid=matrix(0,nrow=nmbrededeltas,ncol=nmbrededeltas)
  
  ##################################################################################################
  # First input
  i = inputs[1]
  
  if (nboot > 0){
    lqidb=array(0,dim=c(nmbrededeltas,nmbrededeltas,nboot)) # pour bootstrap 
    sqhatb=NULL
  }
  
  ## definition of local variables
    Loi.Entree=InputDistributions[[i]]
  
    ############################
    ## definition of local variable vdd
    ############################
    if(!samedelta){
      # In this case, the  mean perturbation will be new_mean = mean+sigma*delta
      # The mean and the standard deviation sigma of the input distribution must be stored in the third place of the list defining the input distribution.
      moy=Loi.Entree[[3]][1]
      sigma=Loi.Entree[[3]][2]
      vdd1=moy+deltasvector*sigma
    } else {
      vdd1=deltasvector
    }
    
    ##########################
    ### The next part does, for each kind of input distribution
    # 	Solve with respect to lambda the following equation :
    #		(1) Mx'(lambda)/Mx(lambda)-delta=0 {See section 3, proposition 3.2}
    #	One can note that (1) is the derivative with respect to lambda of
    #		(2) log(Mx(lambda))-delta*lambda
    # 	Function (2) is concave, therefore its optimisation is theoritically easy.
    #
    #	=> One obtains an unique lambda solution
    #	
    #	Then the density ratio is computed and summed, allowing the estimation of q_i_delta. 
    #	lqid is a vector of same length as  deltas_vector
    #	sigma2_i_tau_N, the estimator of the variance of the estimator of P_i_delta (see lemma 2.1) is also computed
    #
    #	Implemented cases: Gaussian, Uniform, Triangle, Left Trucated Gaussian, Left Truncated Gumbel
    if ( Loi.Entree[[1]] =="norm"||Loi.Entree[[1]] =="lnorm"){
      # if the input is a Gaussian,solution of equation (1) is trivial
      mu1=Loi.Entree[[2]][1]
      sigma1=Loi.Entree[[2]][2]
      phi1=function(tau){mu1*tau+(sigma1^2*tau^2)/2}
      
      vlambda1=(vdd1-mu1)/sigma1^2
    }	# end for Gaussian input, mean twisting
    
    if (  Loi.Entree[[1]] =="unif"){
      # One will not minimise directly log(Mx)-lambda*delta due to numerical problems emerging when considering
      #	exponential of large numbers (if b is large, problems can be expected) 
      # Instead, one note that the optimal lambda for an Uniform(a,b) and a given mean delta is the same that for
      #	an Uniform(0,b-a) and a given mean (delta-a).
      # Function phit corresponding to the log of the Mgf of an U(0,b-a) is implemented; 
      #	the function expm1 is used to avoid numerical troubles when tau is small.
      # Function gt allowing to minimise phit(tau)-(delta-a)*tau is also implemented.
      a1=Loi.Entree[[2]][1]
      b1=Loi.Entree[[2]][2]
      m1=(a1+b1)/2
      
      Mx1=function(tau){
        if (tau==0){ 1 }
        else {(exp(tau*b1)-exp(tau*a1) )/ ( tau * (b1-a1))}
      }
      phi1=function(tau){
        if(tau==0){0}
        else { log ( Mx1(tau))}
      }	
      phit=function(tau){
        if (tau==0){0}
        else {log(expm1(tau*(b1-a1)) / (tau*(b1-a1)))}
      }
      gt=function(tau,delta){ 
        phit(tau) -(delta-a1)*tau
      }
      vlambda1=c();
      for (l in 1:nmbrededeltas){
        tm=nlm(gt,0,vdd1[l])$estimate					
        vlambda1[l]=tm					
      } 	
    }	# end for Uniform input, mean twisting
    
    if (  Loi.Entree[[1]] =="triangle"){
      # One will not minimise directly log(Mx)-lambda*delta due to numerical problems emerging when considering
      #	exponential of large numbers (if b is large, problems can be expected) 
      # Instead, one note that the optimal lambda for an Triangular(a,b,c) and a given mean delta is the same that for
      #	an Uniform(0,b-a,c-a) and a given mean (delta-a).
      # Function phit corresponding to the log of the Mgf of an Tri(0,b-a,c-a) is implemented;
      #	One can note that phit=log(Mx)+lambda*a
      # Function gt allowing to minimise phit(tau)-(delta-a)*tau is also implemented..
      a1=Loi.Entree[[2]][1]
      b1=Loi.Entree[[2]][2]
      c1=Loi.Entree[[2]][3] # reminder: c is between a and b
      m1=(a1+b1+c1)/3	
      
      Mx1=function(tau){
        if (tau !=0){
          dessus=(b1-c1)*exp(a1*tau)-(b1-a1)*exp(c1*tau)+(c1-a1)*exp(b1*tau)
          dessous=(b1-c1)*(b1-a1)*(c1-a1)*tau^2
          return ( 2*dessus/dessous)
        } else {
          return (1)
        }
      }
      phi1=function(tau){return (log (Mx1(tau)))}
      
      phit=function(tau){
        if(tau!=0){
          dessus=(a1-b1)*expm1((c1-a1)*tau)+(c1-a1)*expm1((b1-a1)*tau)
          dessous=(b1-c1)*(b1-a1)*(c1-a1)*tau^2
          return( log (2*dessus/dessous) )
        } else { return (0)}
      }
      gt=function(tau,delta){ 
        phit(tau)-(delta-a1)*tau
      }
      vlambda1=c();
      for (l in 1:nmbrededeltas){
        tm=nlm(gt,0,vdd1[l])$estimate					
        vlambda1[l]=tm					
      } 
    }	# End for Triangle input, mean twisting
    
    if (  Loi.Entree[[1]] =="tnorm"){
      # The case implemented is the left truncated Gaussian
      # Details :	First, the constants defining the distribution (mu, sigma, min) are extracted.
      #		Then, the function g is defined as log(Mx(tau))-delta*tau. 
      #		The function phi has an explicit expression in this case.
      mu1=Loi.Entree[[2]][1]
      sigma1=Loi.Entree[[2]][2]
      min1=Loi.Entree[[2]][3]
      
      phi1=function(tau){
        mpls2=mu1+tau*sigma1^2
        Fa1=pnorm(min1,mu1,sigma1)		
        Fia1=pnorm(min1,mpls2,sigma1)
        lMx=mu1*tau+1/2*sigma1^2*tau^2 - (1-Fa1) + (1-Fia1)
        return(lMx)
      }
      
      g=function(tau,delta){ 
        if (tau == 0 ){	return(0)
        } else {	return(phi1(tau) -delta*tau)}
      }
      vlambda1=c();
      for (l in 1:nmbrededeltas){
        tm=nlm(g,0,vdd1[l])$estimate					
        vlambda1[l]=tm					
      } 
    } 	# End for left truncated Gaussian, mean twisting
    
    
    
    ##################################################################################################
    # Second input
    i = inputs[2]
    
    ## definition of local variables
    Loi.Entree=InputDistributions[[i]]
    
    ############################
    ## definition of local variable vdd
    ############################
    if(!samedelta){
      # In this case, the  mean perturbation will be new_mean = mean+sigma*delta
      # The mean and the standard deviation sigma of the input distribution must be stored in the third place of the list defining the input distribution.
      moy=Loi.Entree[[3]][1]
      sigma=Loi.Entree[[3]][2]
      vdd2=moy+deltasvector*sigma
    } else {
      vdd2=deltasvector
    }
    
    ##########################
    ### The next part does, for each kind of input distribution
    # 	Solve with respect to lambda the following equation :
    #		(1) Mx'(lambda)/Mx(lambda)-delta=0 {See section 3, proposition 3.2}
    #	One can note that (1) is the derivative with respect to lambda of
    #		(2) log(Mx(lambda))-delta*lambda
    # 	Function (2) is concave, therefore its optimisation is theoritically easy.
    #
    #	=> One obtains an unique lambda solution
    #	
    #	Then the density ratio is computed and summed, allowing the estimation of q_i_delta. 
    #	lqid is a vector of same length as  deltas_vector
    #	sigma2_i_tau_N, the estimator of the variance of the estimator of P_i_delta (see lemma 2.1) is also computed
    #
    #	Implemented cases: Gaussian, Uniform, Triangle, Left Trucated Gaussian, Left Truncated Gumbel
    if ( Loi.Entree[[1]] =="norm"||Loi.Entree[[1]] =="lnorm"){
      # if the input is a Gaussian,solution of equation (1) is trivial
      mu2=Loi.Entree[[2]][1]
      sigma2=Loi.Entree[[2]][2]
      phi2=function(tau){mu2*tau+(sigma2^2*tau^2)/2}
      
      vlambda2=(vdd2-mu2)/sigma2^2
    }	# end for Gaussian input, mean twisting
    
    if (  Loi.Entree[[1]] =="unif"){
      # One will not minimise directly log(Mx)-lambda*delta due to numerical problems emerging when considering
      #	exponential of large numbers (if b is large, problems can be expected) 
      # Instead, one note that the optimal lambda for an Uniform(a,b) and a given mean delta is the same that for
      #	an Uniform(0,b-a) and a given mean (delta-a).
      # Function phit corresponding to the log of the Mgf of an U(0,b-a) is implemented; 
      #	the function expm1 is used to avoid numerical troubles when tau is small.
      # Function gt allowing to minimise phit(tau)-(delta-a)*tau is also implemented.
      a2=Loi.Entree[[2]][1]
      b2=Loi.Entree[[2]][2]
      m2=(a2+b2)/2
      
      Mx2=function(tau){
        if (tau==0){ 1 }
        else {(exp(tau*b2)-exp(tau*a2) )/ ( tau * (b2-a2))}
      }
      phi2=function(tau){
        if(tau==0){0}
        else { log ( Mx2(tau))}
      }	
      phit=function(tau){
        if (tau==0){0}
        else {log(expm1(tau*(b2-a2)) / (tau*(b2-a2)))}
      }
      gt=function(tau,delta){ 
        phit(tau) -(delta-a2)*tau
      }
      vlambda2=c();
      for (l in 1:nmbrededeltas){
        tm=nlm(gt,0,vdd2[l])$estimate					
        vlambda2[l]=tm					
      } 	
    }	# end for Uniform input, mean twisting
    
    if (  Loi.Entree[[1]] =="triangle"){
      # One will not minimise directly log(Mx)-lambda*delta due to numerical problems emerging when considering
      #	exponential of large numbers (if b is large, problems can be expected) 
      # Instead, one note that the optimal lambda for an Triangular(a,b,c) and a given mean delta is the same that for
      #	an Uniform(0,b-a,c-a) and a given mean (delta-a).
      # Function phit corresponding to the log of the Mgf of an Tri(0,b-a,c-a) is implemented;
      #	One can note that phit=log(Mx)+lambda*a
      # Function gt allowing to minimise phit(tau)-(delta-a)*tau is also implemented..
      a2=Loi.Entree[[2]][1]
      b2=Loi.Entree[[2]][2]
      c2=Loi.Entree[[2]][3] # reminder: c is between a and b
      m2=(a2+b2+c2)/3	
      
      Mx2=function(tau){
        if (tau !=0){
          dessus=(b2-c2)*exp(a2*tau)-(b2-a2)*exp(c*tau)+(c2-a2)*exp(b2*tau)
          dessous=(b2-c2)*(b2-a2)*(c2-a2)*tau^2
          return ( 2*dessus/dessous)
        } else {
          return (1)
        }
      }
      phi2=function(tau){return (log (Mx2(tau)))}
      
      phit=function(tau){
        if(tau!=0){
          dessus=(a2-b2)*expm1((c2-a2)*tau)+(c2-a2)*expm1((b2-a2)*tau)
          dessous=(b2-c2)*(b2-a2)*(c2-a2)*tau^2
          return( log (2*dessus/dessous) )
        } else { return (0)}
      }
      gt=function(tau,delta){ 
        phit(tau)-(delta-a2)*tau
      }
      vlambda2=c();
      for (l in 1:nmbrededeltas){
        tm=nlm(gt,0,vdd2[l])$estimate					
        vlambda2[l]=tm					
      } 
    }	# End for Triangle input, mean twisting
    
    if (  Loi.Entree[[1]] =="tnorm"){
      # The case implemented is the left truncated Gaussian
      # Details :	First, the constants defining the distribution (mu, sigma, min) are extracted.
      #		Then, the function g is defined as log(Mx(tau))-delta*tau. 
      #		The function phi has an explicit expression in this case.
      mu2=Loi.Entree[[2]][1]
      sigma2=Loi.Entree[[2]][2]
      min2=Loi.Entree[[2]][3]
      
      phi2=function(tau){
        mpls2=mu2+tau*sigma2^2
        Fa2=pnorm(min2,mu2,sigma2)		
        Fia2=pnorm(min2,mpls2,sigma2)
        lMx=mu2*tau+1/2*sigma2^2*tau^2 - (1-Fa2) + (1-Fia2)
        return(lMx)
      }
      
      g=function(tau,delta){ 
        if (tau == 0 ){	return(0)
        } else {	return(phi2(tau) -delta*tau)}
      }
      vlambda2=c();
      for (l in 1:nmbrededeltas){
        tm=nlm(g,0,vdd2[l])$estimate					
        vlambda2[l]=tm					
      } 
    } 	# End for left truncated Gaussian, mean twisting
    
    
    ###########################################################################################
    
    ###############################################################
    ############# Computation of q_i_delta for the mean twisting
    ###############################################################
    
    
    for (K1 in 1:nmbrededeltas){
      for (K2 in 1:nmbrededeltas){
        if ((vdd1[K1]!=0) & (vdd2[K2]!=0)){
          res=NULL ; respts=NULL
          pti1=phi1(vlambda1[K1])
          pti2=phi2(vlambda2[K2])
          for (j in 1:nmbredepoints){	
            res[j]=exp(vlambda1[K1]*xs[j,inputs[1]]-pti1 + vlambda2[K2]*xs[j,inputs[2]]-pti2)
            respts[j]=exp(vlambda1[K1]*x[j,inputs[1]]-pti1 + vlambda2[K2]*x[j,inputs[2]]-pti2)
          }
          sum_res = sum(res)
          kid = 1
          res1 = res[1]
          res2 = res1/sum_res
          while (res2 < order){
            kid = kid + 1
            res1 = res1 + res[kid]
            res2 = res1/sum_res
          }
          if (bias){ lqid[K1,K2] = mean(y[y >= ys$x[kid]]) # ys$x[kid] = quantile
          } else lqid[K1,K2] = mean(y * respts * ( y >= ys$x[kid] ) / (1-order))
        } else lqid[K1,K2] = sqhat
      }
    }
    
    ###############################################################
    ##########              BOOTSTRAP                ###############
    ###############################################################
    
    if (nboot >0){
      for (b in 1:nboot){
        ib <- sample(1:length(y),replace=TRUE)
        xb <- x[ib,]
        yb <- y[ib]
        
        quantilehatb <- quantile(yb,order) # quantile estimate
        sqhatb <- c(sqhatb, mean( yb * (yb >= quantilehatb) / ( 1 - order ))) # superquantile estimate
        
        ysb = sort(yb,index.return=T) # ordered output
        xsb = xb[ysb$ix,] # inputs ordered by increasing output
        
        for (K1 in 1:nmbrededeltas){
          for (K2 in 1:nmbrededeltas){
            if ((vdd1[K1]!=0) & (vdd2[K2]!=0)){
              res=NULL ; respts=NULL
              pti1=phi1(vlambda1[K1])
              pti2=phi2(vlambda2[K2])
              for (j in 1:nmbredepoints){	
                res[j]=exp(vlambda1[K1]*xsb[j,inputs[1]]-pti1 + vlambda2[K2]*xsb[j,inputs[2]]-pti2)
                respts[j]=exp(vlambda1[K1]*xb[j,inputs[1]]-pti1 + vlambda2[K2]*xb[j,inputs[2]]-pti2)
              }
              sum_res = sum(res)
              kid = 1
              res1 = res[1]
              res2 = res1/sum_res
              while (res2 < order){
                kid = kid + 1
                res1 = res1 + res[kid]
                res2 = res1/sum_res
              }
              if (bias){ lqidb[K1,K2,b] = mean(yb[yb >= ysb$x[kid]]) # ysb$x[kid] = quantile
              } else lqidb[K1,K2,b] = mean(yb * respts * (yb >= ysb$x[kid] ) / (1-order))
            } else lqidb[K1,K2,b] = sqhatb[b]
          }
        }
      } # end of bootstrap loop
    }
    
  ######################################################################################
  
  ########################################
  ## Plugging estimator of the S_i_delta indices
  ########################################	
    
  for (i in 1:nmbrededeltas){
    for (j in 1:nmbrededeltas){
      J[i,j]=lqid[i,j]
      if (percentage==FALSE) I[i,j]=transinverse(lqid[i,j],sqhat,sqhat)
      else I[i,j]=lqid[i,j]/sqhat-1
      
      if (nboot > 0){
        # ICI = PLI bootstrap including or excluding sampling uncertainty
        # JCI = quantile bootstrap 
        sqinf <- quantile(lqidb[i,j,],(1-conf)/2)
        sqsup <- quantile(lqidb[i,j,],(1+conf)/2)
        JCIinf[i,j]=sqinf
        JCIsup[i,j]=sqsup
        sqb <- mean(sqhatb)
        if (percentage==FALSE){
          if (bootsample){
            ICIinf[i,j]=transinverse(sqinf,sqb,sqb)
            ICIsup[i,j]=transinverse(sqsup,sqb,sqhat)
          } else{
            ICIinf[i,j]=quantile(transinverse(lqidb[i,j,],sqhatb,sqhatb),(1-conf)/2)
            ICIsup[i,j]=quantile(transinverse(lqidb[i,j,],sqhatb,sqhatb),(1+conf)/2)
          }
        } else {
          if (bootsample){
            ICIinf[i,j]=sqinf/sqb-1
            ICIsup[i,j]=sqsup/sqb-1
          } else{
            ICIinf[i,j]=quantile((lqidb[i,j,]/sqhatb-1),(1-conf)/2)
            ICIsup[i,j]=quantile((lqidb[i,j,]/sqhatb-1),(1+conf)/2)
          }
        }
      }
    }
  }

  res <- list(PLI = I, PLICIinf = ICIinf, PLICIsup = ICIsup, quantile = J, quantileCIinf = JCIinf, quantileCIsup = JCIsup)
  
  
  return(res)
}
