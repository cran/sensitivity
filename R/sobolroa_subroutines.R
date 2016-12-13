# -------------------------------------------------------------------
#
# SOBOLROAUC/SOBOLROALHS SUBROUTINES
# 
# -------------------------------------------------------------------

# Purpose:
# 
#   HADAM calculates the Hadamard matrix of the Galois field GF(q)
# 
# Modified:
# 
#   December 9, 2016
# 
# Author:
# 
#   Laurent Gilquin
# 
# Parameters:
# 
#   Input, int Q, a prime number 
#
#   Output, matrix HADAM, the multiplication table of GF(q)

sobolroa.Hadam=function(q){
  
  S <- seq(0,(q-1))
  
  return(S%*%t(S)%%q)
}

# -----------------------------------------------------------------------
# Purpose:
# 
#   PERMUTATIONS generates all permutations of length N.
# 
# Modified:
# 
#   December 9, 2016
# 
# Author:
# 
# Credit goes to "Museful" from stackoverflow.
# 
# Parameters:
# 
#   Input, int N, a positive integer
#
#   Output, matrix PERMUTATIONS, the matrix storing all the permutations (N! rows).

sobolroauc.permutations <- function(n){
  if(n==1){
    return(matrix(1))
  } else {
    sp <- sobolroauc.permutations(n-1)
    p <- nrow(sp)
    A <- matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}


# ---------------------------------------------------------------------
# Purpose:
# 
#   SIMPLEX_ORD returns the vertices of the unit ordered simplex in the d-hypercube
# 
# Modified:
# 
#   December 9, 2016
# 
# Author:
# 
#   Laurent Gilquin
# 
# Parameters:
# 
#   Input, int D, the dimension of the hypercube
#
#   Output, matrix SIMPLEX_ORD, a matrix storing by row the vertices sorted in
#   ascending dimensions  

sobolroauc.simplex_ord=function(d){
  
  M <- matrix(NA,nrow=d+1,ncol=d)
  for (i in 1:d){
    M[,i] <- c(rep(0,d+1-i),rep(1,i))
  }
  return(M)
}


# ---------------------------------------------------------------------------
# Purpose:
# 
#   SIMPLEX_ID returns the vertices ranks of the d! simplices contained in the d-hypercube.
# 
# Modified:
# 
#   December 9, 2016
# 
# Author:
# 
#   Laurent Gilquin
# 
# Parameters:
# 
#   Input, int D, the dimension of the hypercube
#
#   Output, matrix SIMPLEX_ID, the matrix storing by row the vertices ranks (d! rows) 

sobolroauc.simplex_id=function(d){
  
  ref <- sobolroauc.simplex_ord(d)
  index <- sobolroauc.permutations(d)
  x <- c(0,1)
  X <- expand.grid(as.data.frame(matrix(rep(x,d),ncol=d)))
  X <- data.frame(t(unname(X)))
  M <- data.frame(t(matrix(ref[,c(index)],ncol=d)))
  vectors <- matrix(match(M,X),ncol=d+1,byrow=TRUE)
  return(vectors)
}

# --------------------------------------------------------------------------
# Purpose:
# 
#   ELEMENT subdivises a d-hypercube in m^d sub-hypercubes of same volume.
# 
# Modified:
# 
#   December 9, 2016
# 
# Author:
# 
#   Laurent Gilquin
# 
# Parameters:
# 
#   Input, int D, the dimension of the hypercube
#
#   Input, int M, the length of the margins subdivision 
#
#   Output, matrix ELEMENT, matrix storing by row the vertices ranks of each sub-hypercube 


sobolroauc.element=function(m,d){
  
  l <- seq(1,(m-1))
  ind <- seq(0,(m-2))
  for (i in 1:(d-1)){
    l <- c(outer(l,ind*m^i,"+"))
  }
  for (i in 0:(d-1)){
    l <- c(outer(l,c(0,m^i),"+"))
  }
  M <- matrix(l,ncol=2^d)
  return(M)
}

# -----------------------------------------------------------------------
# Purpose:
# 
#   SIMPLEXATE subdivises one or multiple sub-hypercubes in simplices
#
# Modified:
# 
#   December 9, 2016
# 
# Author:
# 
#   Laurent Gilquin
# 
# Parameters:
# 
#   Input, matrix V, the vertices ranks of the simplices
#
#   Input, matrix ELE, the vertices ranks of each sub-hypercube
#
#   Output, matrix SIMPLEXATE, matrix storing by row the vertices ranks of each simplices


sobolroauc.simplexate=function(v,ele){
  
  vector <- as.vector(t(v))
  simplex <- matrix(t(ele[,c(vector)]),ncol=ncol(v),byrow=TRUE)
  return(simplex)
}

# ------------------------------------------------------------------------------
# Purpose:
# 
#   SIMPLEX_TRUE returns among a set the simplices those satisfaying the full set of linear ordered constraints
#
# Modified:
# 
#   December 9, 2016
# 
# Author:
# 
#   Laurent Gilquin
# 
# Parameters:
# 
#   Input, matrix SIMPLEX, vertices ranks of the simplices
#
#   Input, matrix NODES, the grid of points
#
#   Output, matrix SIMPLEX_TRUE, matrix storing by row the vertices ranks of each simplex
#   satisfaying the full set of linear ordered constraints


sobolroauc.simplex_true=function(simplex,nodes){
  
  d <- ncol(nodes)
  vec <- rowSums(matrix(c(nodes[,-d])<=c(nodes[,-1]),ncol=d-1))
  ind <- which(vec<(d-1)) 
  M <- matrix(c(simplex)%in%ind,ncol=ncol(simplex))
  keep <- which(rowSums(M)==0)
  return(simplex[keep,]) 
}

# ---------------------------------------------------------------------------
# Purpose:
# 
#   SIMPLEX_CREATE subdivises the unit ordered simplex of a d-hypercube into simplices
#
# Modified:
# 
#   December 9, 2016
# 
# Author:
# 
#   Laurent Gilquin
# 
# Parameters:
# 
#   Input, int M, the length of the margins subdivision
#
#   Input, int D, the dimension of the d-hypercube
#
#   Output, list SIMPLEX_CREATE, list containing:
#   matrix S, the matrix storing by row the vertices ranks of the simplices
#   matrix P, the matrix containing the points created with the subdivision 

sobolroauc.simplex_create=function(m,d){
  
  x <- seq(0,1,len=m)
  nodes <- expand.grid(as.data.frame(matrix(rep(x,d),ncol=d)))
  nodes <- as.matrix(unname(nodes))
  vectors <- sobolroauc.simplex_id(d)
  ele <- sobolroauc.element(m,d)
  simplex <- sobolroauc.simplexate(vectors,ele)
  simplex <- sobolroauc.simplex_true(simplex,nodes)
  return(list(s=simplex,p=nodes))
}
