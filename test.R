require(sensitivity)
   n <- 1000
     X1 <- data.frame(matrix(runif(3 * n), nrow = n))
     X2 <- data.frame(matrix(runif(3 * n), nrow = n))
     
     # sensitivity analysis
     x=sobolCert(model=function(X) { list(out=X[1]+2*X[2]+X[3]+.001*runif(1),err=.01); }, 
                 X1, X2, conf=.99, lambda0=.1, h=.1, nboot=30)
     print(x)
     
     x=sobolCert(model=NULL, X1=NULL, X2=NULL, conf=.95)
     print(x)

