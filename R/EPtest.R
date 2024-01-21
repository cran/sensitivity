EPtest = function(X, y, u = NULL, doe = NULL, Kdoe = 10, tau = 0.1){

  ###################################
  ## Description of input parameters
  ###################################
  	
  # X is a matrix or data.frame that contains the numerical inputs as columns.
  
  # y is the output vector (numeric).
  
  # u is the vector of indices of the columns of X for which we want to test the significance. 
  #     If Su denotes the Sobol index associated to the input columns i in u, and S the total Sobol index,
  #     the null hypothesis is H0 : Su = S.
  
  # doe is the design of experiment on which the empirical process xi is evaluated. 
  
  # If doe is null and Kdoe is specified, the design of experiment is taken as Kdoe points drawn uniformly 
  #     on the range of the input matrix X.
  
  # tau in (0;1) is a regularization parameter to approximate the limit chi2 distribution of the 
  #     test statistics under H0. Once the asymptotic covariance matrix Sigma of the process xi 
  #     is estimated, only the eigenvalues of Sigma above tau times its spectral radius are kept. 
  
  ###################################
  ## Description of output parameters
  ###################################

  # statistics is the test statistics expected to follow a chi-squared distribution under the null hypothesis.
  
  # ddl is the number of degrees of freedom used for the test.
  
  # pv is the test p-value
  
  ########################################
  ## Local variables	
  ########################################
  
  lu = length(u);
  if(is.vector(X)){
    X = matrix(X, ncol = 1)
  }
  X = as.matrix(X);
  n = nrow(X);
  nc = ncol(X);
  ubar = setdiff(1:nc, u)
  
  ########################################
  ## vector comparison function
  ########################################
  
  leq = function(Z, z){
    Z = as.matrix(Z)
    if(length(z) != ncol(Z)){
      stop("dimension issue")
    }
    return(apply(Z, 1, function(v){all(v <= z)}))
  }
  
  ########################################
  ## empirical moments function
  ########################################
  
  m = function(x, k, u = 1:length(x)){
    if(length(u) == 0){
      return(mean(y**k))
    }
    return(mean(y**k * leq(X[, u], x[u])))
  }
  
  ########################################
  ## eta process function
  ########################################
  
  eta = function(x){
    return(c(m(x, k = 1), 
             m(x, k = 1, u = u), 
             m(x, k = 0, u = ubar)))
  }
  
  ########################################
  ## xi process function
  ########################################
  
  xi = function(xx){
    nr = nrow(xx)
    xi = rep(NA, nr)
    for(i in 1:nr){
      e = eta(xx[i, ])
      xi[i] = e[1] - e[2]*e[3]
    }
    return(xi)
  }
  
  ########################################
  ## Omega matrix function
  ########################################
  
  Omega = function(x, xp){
    p = pmin(x, xp)
    pu = pmin(x[u], xp[u])
    M = matrix(NA, nrow = 3, ncol = 3)
    M[1, 1] = m(p, k = 2)
    xx = x ; xx[u] = pu ; M[1, 2] = m(xx, k = 2)
    xx = p ; xx[u] = x[u] ; M[1, 3] = m(xx, k = 1)
    xx = xp ; xx[u] = pu ; M[2, 1] = m(xx, k = 2)
    M[2, 2] = m(p, k = 2, u = u)
    xx = xp ; xx[u] = x[u] ; M[2, 3] = m(xx, k = 1)
    xx = p ; xx[u] = xp[u] ; M[3, 1] = m(xx, k = 1)
    xx = x ; xx[u] = xp[u] ; M[3, 2] = m(xx, k = 1)
    M[3, 3] = m(p, k = 0, u = ubar)
    return(M - eta(x) %*% t(eta(xp)))
  }
  
  ########################################
  ## Sigma matrix function
  ########################################
  
  sigma = function(x, xp){
    v = c(1, -m(x, k = 0, u = ubar), -m(x, k = 1, u = u))
    vp = c(1, -m(xp, k = 0, u = ubar), -m(xp, k = 1, u = u))
    return(t(v) %*% Omega(x, xp) %*% vp)
  }
  
  ########################################
  ## Generate the design of experiment on which the process xi is to be evaluated
  ########################################
  
  if(is.null(doe)){
    K = Kdoe
    if(nc == 1){
      minX = min(X)
      maxX = max(X)
      xx = rep(minX, K) + (maxX - minX)*runif(K)
      xx = matrix(xx, ncol = 1)
    }
    if(nc >1){
      RU = matrix(runif(nc*K), ncol = nc)
      minX = apply(X, 2, min)
      maxX = apply(X, 2, max)
      xx = matrix(rep(1, K), ncol = 1) %*% matrix(minX, nrow = 1) + RU %*% diag(maxX - minX)
    }
  }
  
  if(!is.null(doe)){
    xx = as.matrix(doe)
    K = nrow(xx)
  }

  ########################################
  ## Empirical estimator of the asymptotic covariance matrix of xi
  ########################################
  
  Sigma = matrix(NA, K, K)
  for(i in 1:K){
    for(j in i:K){
      Sigma[i, j] = sigma(xx[i, ], xx[j, ])
      if(j != i){
        Sigma[j, i] = Sigma[i, j]
      }
    }
  }
  
  ########################################
  ## regularization step
  ########################################
  
  e = eigen(Sigma, symmetric = TRUE)
  lambda = e$values
  nl = sum(lambda > tau * n^(-1/3) * lambda[1])
  lambda = lambda[1:nl]
  
  ########################################
  ## test statistics and p-value
  ########################################
  
  Q = n * sum((t(xi(xx)) %*% as.matrix(e$vectors[, 1:nl], ncol = nl))**2 / lambda)
  pv = 1 - pchisq(Q, df = nl)
  
  return(list("statistics" = Q, "ddl" = nl, "p-value" = pv))
}