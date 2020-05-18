discrepancyCriteria_cplus <- function(design, type = "all") {
  
  X <- as.matrix(design)
  dimension <- dim(X)[2]
  n <- dim(X)[1]
  if (n < dimension) {
    stop("Warning : the number of points is lower than the dimension.")
  }
  if (min(X) < 0 || max(X) > 1) {
    warning("The design is rescaling into the unit cube [0,1]^d.")
    M <- apply(X, 2, max)
    m <- apply(X, 2, min)
    for (j in 1:dim(X)[2]) {
      X[, j] <- (X[, j] - m[j])/(M[j] - m[j])
    }
  }
  R <- list()
  DisC2 <- FALSE
  DisL2 <- FALSE
  DisL2star <- FALSE
  DisM2 <- FALSE
  DisS2 <- FALSE
  DisW2 <- FALSE
  if (length(type) == 1 && type == "all") {
    type <- c("C2", "L2", "L2star", "M2", "S2", "W2")
  }
  for (i in 1:length(type)) {
    type_ <- type[i]
    switch(type_, C2 = {
      DisC2 <- TRUE
    }, L2 = {
      DisL2 <- TRUE
    }, L2star = {
      DisL2star <- TRUE
    }, M2 = {
      DisM2 <- TRUE
    }, S2 = {
      DisS2 <- TRUE
    }, W2 = {
      DisW2 <- TRUE
    })
  }
  if (DisC2 == TRUE) {
    
    P <- 1 + 0.5 * abs(X - 0.5) - 0.5 * (abs(X - 0.5)^2)
    s1 <- DisC2_Rowprod(t(P),dimension)
    s2 <- DisC2_Crossprod(c(t(X)),dimension)
    R <- c(R, DisC2 = sqrt(((13/12)^dimension) - ((2/n) * s1) + ((1/n^2) * s2)))
  }
  
  if (DisL2 == TRUE) {
    
    P <- X*(1-X)
    s1 <- DisL2_Rowprod(t(P),dimension)
    s2 <- DisL2_Crossprod(c(t(X)),dimension)
    R <- c(R, DisL2 = sqrt(12^(-dimension) - (((2^(1 - dimension))/n) * s1) + ((1/n^2) * s2)))        
  }
  
  if (DisL2star == TRUE) {
    
    dL2 <- DisL2star_Crossprod(t(X),dimension)
    R <- c(R, DisL2star = sqrt(3^(-dimension) + dL2))
  }
  
  if (DisM2 == TRUE) {
    
    P <- 3-X^2
    s1 <- DisM2_Rowprod(t(P),dimension)
    s2 <- DisM2_Crossprod(c(t(X)),dimension)
    R <- c(R, DisM2 = sqrt(((4/3)^dimension) - (((2^(1 - dimension))/n) * s1) + ((1/n^2) * s2)))
    
  }
  
  if (DisS2 == TRUE) {
    
    P <- 1+2*X-2*X^2
    s1 <- DisS2_Rowprod(t(P),dimension)
    s2 <- DisS2_Crossprod(c(t(X)),dimension)
    R <- c(R, DisS2 = sqrt(((4/3)^dimension) - ((2/n) * s1) + ((2^dimension/n^2) * s2)))
  }
  
  if (DisW2 == TRUE) {
    
    s1 <- DisW2_Crossprod(t(X),dimension)
    R <- c(R, DisW2 = sqrt(-(4/3)^dimension + (1/n^2) * s1))
  }
  
  return(R)
}