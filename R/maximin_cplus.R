maximin_cplus <- function(design){
  
  X <- as.matrix(design)
  n <- dim(X)[1]
  # To check the experimental region
  if ( min(X)<0 || max(X)>1 ){
    warning("The design is rescaling into the unit cube [0,1]^d.")
    M <- apply(X,2,max)
    m <- apply(X,2,min)
    for (j in 1:dim(X)[2]){
      X[,j] <- (X[,j]-m[j])/(M[j]-m[j])
    }	
  }
  
  val <- maximin_cpp(X)
  return(val)
}
