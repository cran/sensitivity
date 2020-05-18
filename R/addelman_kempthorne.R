# Addelman and Kempthorne construction of orthogonal arrays
# Heydayat et al. (1999), Orthogonal Arrays : Theory and Applications, p.127-128

# function testing if a matrix is a difference scheme
is_diff_scheme <- function(M,sub_table,q){
  M <- M[,-1]
  d <- ncol(M)
  test <- TRUE
  i <- 1
  j <- 1
  l <- 1
  while (test & (l < (d*(d-1)/2))) {
    i <- i+(j==d)
    j <- ((j==d)*i + (j<d)*j)+1
    test <- length(unique(sub_table[cbind(M[,i],M[,j])+1]))==q
    l <- l+1
  }
  return(test)
}

# outliers test function
test_outlier <- function(OA1,OA2){
  d <- ncol(OA1)
  i <- 1
  j <- 1
  for (l in 1:(d*(d-1)/2)){
    i <- i+(j==d)
    j <- ((j==d)*i + (j<d)*j)+1
    Ind1 <- do.call(base::order,as.data.frame(OA1[,c(i,j)]))
    Ind2 <- do.call(base::order,as.data.frame(OA2[,c(i,j)]))
    test <- max(abs(cor(OA1[Ind1,-c(i,j)],OA2[Ind2,-c(i,j)])))>10^(-1)
    if(test){
      break 
    }
  }
  return(test)
}

# function returning multiplication and addition tables of GF(q)
field_tables <- function(q){
  prime_fact <- primeFactors(q)
  if(length(prime_fact)==1){
    gf_mul <- (matrix(0:(q-1),ncol=1)%*%(0:(q-1)))%%q
    gf_add <- t(t(replicate(q,0:(q-1)))+0:(q-1))%%q
  } else {
    gf_mul <- get(paste("gf",q,"_mul",sep=""))
    gf_add <- get(paste("gf",q,"_add",sep=""))
  }
  return(list(mul=gf_mul,add=gf_add))
}

# function returning the list of coefficients for the Addelman and Kempthorne construction 
field_coeffs <- function(q){
  prime_fact <- primeFactors(q)
  if(length(prime_fact)==1){
    gf_mul <- field_tables(q)$mul
    gf_vu <- primroot(q)
    gf_beta <- which(gf_mul[3,]==1)-1
    gf_gamma <- ((gf_vu-1)*(which(gf_mul[(2*gf_vu)%%q+1,]==1)-1))%%q
    gf_delta <- (gf_vu*gf_beta)%%q
    gf_epsilon <- ((gf_vu-1)*gf_beta)%%q
    gf_coeff_prime <- c(gf_vu,gf_beta,gf_gamma,gf_delta,gf_epsilon)
  } else {
    gf_coeff_prime <- get(paste("gf",q,"_coeff_prime",sep=""))[1,]
  }
  return(gf_coeff_prime)
}

# Addelman and Kempthorne construction
addelman_const <- function(dimension, levels, choice="U"){
  
  d <- dimension
  q <- levels
  
  if(d>q){
    stop("The number of columns of the orthogonal array must be lower or equal to 
         its number of levels.")
  }
  
  coeff <- field_coeffs(q)
  tables <- field_tables(q)
  gf_mul <- tables$mul
  gf_add <- tables$add
  if(choice=="U"){
    diff_scheme <- gf_mul
  } 
  if(choice=="V"){
    diff_scheme <- matrix(gf_add[cbind(c(gf_mul),rep(gf_mul[coeff[2]+1,diag(gf_mul)+1],each=q))+1],nrow=q)
  }
  if(choice=="W"){
    diff_scheme <- matrix(gf_add[cbind(c(gf_mul),rep(gf_mul[coeff[3]+1,diag(gf_mul)+1],q))+1],nrow=q)
  }
  if(choice=="X"){
    X_bis <- matrix(gf_add[cbind(gf_mul[coeff[1]+1,gf_mul+1],rep(gf_mul[coeff[4]+1,diag(gf_mul)+1],each=q))+1],nrow=q)
    diff_scheme <- matrix(gf_add[cbind(c(X_bis),rep(gf_mul[coeff[5]+1,diag(gf_mul)+1],q))+1],nrow=q)
  }
  OA <- matrix(c(gf_add[t(diff_scheme[,sample(q,d)]+1),]),byrow=TRUE,ncol=d)+1
  return(OA)
}
