
shapleyPermEx <- function(model = NULL, Xall, Xset, d, Nv, No, Ni = 3, colnames = NULL, ...) {
  #################################################################################
  # Authors: 
  # First version by: Eunhye Song, Barry L. Nelson, Jeremy Staum (Northwestern University), 2015
  # Modified for 'sensitivity' package by: Bertrand Iooss (EDF R&D), 2017
  #################################################################################
  
  #################################################################################
  # This function implements Algorithm 1 to calculate the Shapley effects
  # by examining all permutations of inputs.
  # It also estimates First order and total Sobol' indices 
  #
  # List of inputs to this function:
  # model: a function, or a model with a predict method, defining the model to analyze
  # Xall(n): a function to generate a n-sample of a d-dimensional input vector
  # Xset(n, Sj, Sjc, xjc): a function to generate a n- sample an input vector corresponding 
  #   to the indices in Sj conditional on the input values xjc with the index set Sjc
  # d: number of inputs
  # Nv: Monte Carlo (MC) sample size to estimate the output variance
  # No: output MC sample size to estimate the cost function
  # Ni: inner MC sample size to estimate the cost function
  #
  # This function requires R package "gtools" 
  #
  #################################################################################
  
  # Generate all d! permutations
  perms <- gtools::permutations(d, d, 1:d)
  m <- nrow(perms)
  
  X <- matrix(NA, ncol=d, nrow=Nv+m*(d-1)*No*Ni)
  if (is.null(colnames)) for (i in 1:d) colnames <- c(colnames,paste("X",i,sep=""))
  colnames(X) <- colnames
  
  X[1:Nv,] <- Xall(Nv)
  
  for (p in 1:m)
  {
    pi <- perms[p,]
    pi_s <- sort(pi,index.return=TRUE)$ix
    
    for (j in (1:(d-1)))
    {
      Sj <- pi[c(1:j)] # set of the 1st-jth elements in pi 
      Sjc <- pi[-c(1:j)] # set of the (j+1)th-dth elements in pi
        
      xjcM <- matrix(Xset(No, Sjc, NULL, NULL),nrow=No) # sampled values of the inputs in Sjc
      for (l in 1:No)
      {
        xjc <- xjcM[l,]
          
        # sample values of inputs in Sj conditional on xjc
        xj <- Xset(Ni, Sj, Sjc, xjc) 
        xx <- cbind(xj, matrix(xjc,nrow=Ni,ncol=length(xjc),byrow=T))
        X[(Nv+(p-1)*(d-1)*No*Ni+(j-1)*No*Ni+(l-1)*Ni+1):(Nv+(p-1)*(d-1)*No*Ni+(j-1)*No*Ni+l*Ni),] <- xx[,pi_s]
      }
    }
  }
  
  sh <- list(model = model, Xall = Xall, xset = Xset, d = d, Nv = Nv, No = No, Ni = Ni, X = X, 
             colnames = colnames, perms = perms, call = match.call())
  class(sh) <- "shapleyPermEx"
    
  if (!is.null(sh$model)) {
    response(sh, ...)
    tell(sh)
  }
    
  return(sh)
}


tell.shapleyPermEx <- function(x, y = NULL, return.var = NULL, ...) {
  id <- deparse(substitute(x))
  
  if (! is.null(y)) {
    x$y <- y
  } else{ 
    y <- x$y
    if (is.null(x$y)) stop("y not found")
  }
  
  d <- x$d
  
  # Initialize Shapley value for all players
  Sh <- rep(0, d) ; Sh2 <- rep(0, d)
  
  # Initialize main and total (Sobol) effects for all players
  Vsob <- rep(0, d) ; Vsob2 <- rep(0, d)
  Tsob <- rep(0, d) ; Tsob2 <- rep(0, d)
  
  # Estimate Var[Y] 
  Y <- y[1:x$Nv] ; y <- y[-(1:x$Nv)]
  EY <- mean(Y)
  VarY <- var(Y)
  
  # Estimate Shapley effects
  m <- nrow(x$perms)
  for (p in 1:m)
  {
    pi <- x$perms[p,]
    prevC <- 0
    for (j in 1:d)
    {
      if (j == d)
      {    
        Chat <- VarY
        del <- Chat - prevC
        Vsob[pi[j]] <- Vsob[pi[j]] + prevC # first order effect
        Vsob2[pi[j]] <- Vsob2[pi[j]] + prevC^2
      }
      else
      {
        cVar <- NULL
        for (l in 1:x$No)
        {
          Y <- y[1:x$Ni] ; y <- y[-(1:x$Ni)]
          cVar <- c(cVar,var(Y))
        }
        Chat <- mean(cVar)
        del <- Chat - prevC
      }
      
      Sh[pi[j]] <- Sh[pi[j]] + del
      Sh2[pi[j]] <- Sh2[pi[j]] + del^2
      
      prevC <- Chat
      
      if (j == 1){
        Tsob[pi[j]] <- Tsob[pi[j]] + Chat # Total effect
        Tsob2[pi[j]] <- Tsob2[pi[j]] + Chat^2
      }
    }
  }
  Sh <- Sh / m / VarY
  Sh2 <- Sh2 / m / VarY^2
  ShSE <- sqrt((Sh2 - Sh^2) / m)
  
  Vsob <- Vsob / (m/d) / VarY # averaging by number of permutations with j=d-1
  Vsob2 <- Vsob2 / (m/d) / VarY^2
  VsobSE <- sqrt((Vsob2 - Vsob^2) / (m/d))
  Vsob <- 1 - Vsob 
  Vsob2 <- 1 - Vsob2 
  
  Tsob <- Tsob / (m/d) / VarY # averaging by number of permutations with j=1
  Tsob2 <- Tsob2 / (m/d) / VarY^2
  TsobSE <- sqrt((Tsob2 - Tsob^2) / (m/d))
  
  Shapley <- data.frame(Sh,row.names=x$colnames)
  names(Shapley) <- c("original")
  Vsobol <- data.frame(Vsob,row.names=x$colnames)
  names(Vsobol) <- c("original")
  Tsobol <- data.frame(Tsob,row.names=x$colnames)
  names(Tsobol) <- c("original")
  # Following lines are commented because standard deviation estimates are incorrect
  #Shapley <- data.frame(cbind(Sh,ShSE,Sh-2*ShSE,Sh+2*ShSE),row.names=x$colnames)
  #names(Shapley) <- c("original","std. error", "min. c.i.", "max. c.i.")
  #Vsobol <- data.frame(cbind(Vsob,VsobSE,Vsob-2*VsobSE,Vsob+2*VsobSE),row.names=x$colnames)
  #names(Vsobol) <- c("original","std. error", "min. c.i.", "max. c.i.")
  #Tsobol <- data.frame(cbind(Tsob,TsobSE,Tsob-2*TsobSE,Tsob+2*TsobSE),row.names=x$colnames)
  #names(Tsobol) <- c("original","std. error", "min. c.i.", "max. c.i.")
  
  x$Shapley <- Shapley
  x$SobolS <- Vsobol
  x$SobolT <- Tsobol
  x$V <- VarY
  x$E <- EY
  
  for (i in return.var) {
    x[[i]] <- get(i)
  }
  
  assign(id, x, parent.frame())
}


print.shapleyPermEx <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (!is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    cat("\nShapley' effects:\n")
    print(x$Shapley)
    cat("\nFirst order Sobol' indices:\n")
    print(x$SobolS)
    cat("\nTotal Sobol' indices:\n")
    print(x$SobolT)
  }
}


plot.shapleyPermEx <- function(x, ylim = c(0, 1), ...) {
  if (!is.null(x$y)) {
    pch = c(21, 24, 25)
    nodeplot(x$Shapley, xlim = c(1, x$d + 1), ylim = ylim, pch = pch[1])
    nodeplot(x$SobolS, xlim = c(1, x$d + 1), ylim = ylim, labels = FALSE,
             pch = pch[2], at = (1:x$d)+.2, add = TRUE)
    nodeplot(x$SobolT, xlim = c(1, x$d + 1), ylim = ylim, labels = FALSE,
             pch = pch[3], at = (1:x$d)+.4, add = TRUE)
    legend(x = "topright", legend = c("Shapley effect","First Sobol' index", "Total Sobol' index"), pch = pch)
  }
}

