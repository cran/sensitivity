meanknn <- function(Y,id.sort,n.knn){
  # Estimation of E(Var(Y|Xu))/Var(Y)
  q <- ncol(Y)
  res <- apply(Y,2,function(y){
    mean(apply(matrix(y[id.sort],ncol=n.knn),1,var))/var(y)
  })
  return(unname(res))
}

rescale_inputs <- function(X){
  # Mahalanobis whitening, new X has unit diagonal covariance matrix
  X <- whitening::whiten(X, method="ZCA-cor")
  # We then apply a copula transform
  X <- apply(X,2,rank)/nrow(X)
  return(X)
}

sobolshap_knn <- function(model = NULL, X, id.cat=NULL, U=NULL, method="knn", n.knn=2,
                          return.shap=FALSE, randperm=FALSE, n.perm=1e4, rescale=FALSE, n.limit=2000, noise=FALSE, ...) {
  
  id.cat.detected <- which(sapply(X,class)=="factor")
  if (length(id.cat.detected)){
    # Some inputs are factors
    if (is.null(id.cat)){
      warning("Some inputs are factors but id.cat=NULL, replacing id.cat with the indices of the detected factor inputs")
      id.cat <- id.cat.detected
    }else{
      if (!all(id.cat=id.cat.detected)){
        warning("id.cat is not consistent with detected factors, replacing id.cat with the indices of the detected factor inputs")
        id.cat <- id.cat.detected
      }
    }
  }
  
  if (!is.null(U) & return.shap){
    # return.shap can only be true if U is null (i.e. we generate all combinations)
    warning("We can only compute Shapley values if U is null, switching to return.shap=FALSE")
    return.shap=FALSE
  }
  
  if (is.null(U) & !return.shap){
    warning("U is null, switching to return.shap=TRUE")
    return.shap=TRUE
  }
  
  p <- ncol(X)
  if (is.null(colnames(X))){
    colnames(X) <- paste0("X",1:p)
  }
  
  if (!is.data.frame(X)){
    X <- as.data.frame(X)
  }
  
  x <- list(model = model, X = X, id.cat=id.cat, U=U, method=method, n.knn=n.knn,
            return.shap=return.shap, randperm=randperm, n.perm=n.perm, rescale=rescale, n.limit=n.limit, noise=noise, call = match.call()) 
  class(x) <- "sobolshap_knn"
  
  #calcul of the response for explicit model
  if (! is.null(x$model)) {
    response(x, ...)
    x=tell(x, ...)
  }  
  return(x)
}

estim.sobolshap_knn <- function(data, i=1:nrow(data), q, id.cat,U,method,n.knn,
                                return.shap,randperm,n.perm,rescale,n.limit,noise){
  
  ptot <- ncol(data)
  p <- ptot - q
  X <- data[i,1:p]
  Y <- data[i,(p+1):(ptot),drop=FALSE]
  
  n.cat <- length(id.cat) # number of categorical input variables
  n <- nrow(X) # sample size
  
  # Normalization and encoding of categorical variables
  # ---------------------------------------------------
  if (n.cat){
    id.var <- vector("list",p)
    # Re-arrange X matrix with continuous variables first and rescale if necessary
    id.cont <- setdiff(1:p,id.cat)
    Xtemp <- as.matrix(X[,id.cont,drop=FALSE])
    if (rescale & ncol(Xtemp) > 1){
      Xtemp <- rescale_inputs(Xtemp)
    }
    id.var[id.cont] <- 1:length(id.cont)
    # One-hot encoding of categorical variables, with normalization
    for (i in 1:n.cat){
      df.temp <- X[,id.cat[i],drop=FALSE]; colnames(df.temp) <- "X"
      xx.temp <- model.matrix(~X-1,df.temp)*sqrt(2)/2 # normalization ensures that two different rows will have unit euclidean distance
      colnames(xx.temp) <- paste0("Xf",i,1:ncol(xx.temp))
      id.var[[id.cat[i]]] <- (ncol(Xtemp)+1):(ncol(Xtemp)+ncol(xx.temp))
      Xtemp <- cbind(Xtemp,xx.temp)
    }
    X <- Xtemp
  }else{
    id.var <- as.list(1:p)
    X <- as.matrix(X)
    if (rescale  & ncol(X) > 1){
      X <- rescale_inputs(X)
    }
  }
  
  # Define the U matrix containing the subsets for which
  # we want to compute the Sobol index
  # ------------------------------------------------------
  # order will be NULL unless the user asks for Sobol first-order (order=1) or total (order=0) indices 
  if (is.null(U)){
    order <- NULL
    if (!randperm){
      # We generate all possible subsets
      U <- t(sapply(1:(2^p-2),function(x){as.numeric(intToBits(x)[1:p])}))
      U <- rbind(U,rep(1,p))
    }else{
      # We select the subset with random permutations
      U <- NULL
      perms <- matrix(NA,n.perm,p)
      for (i in 1:n.perm){
        perm <- sample.int(p,p)
        Utemp <- t(sapply(1:p,function(id){
          row <- rep(0,p)
          row[perm[1:id]] <- 1
          return(row)
        }))
        U <- rbind(U,Utemp)
        perms[i,] <- perm
      }
      U <- unique(U) # We remove redundant subsets, if any
    }
  }else{
    if (is.matrix(U)){
      # Only subsets of variables given by the user in U
      order <- NULL
    }else if (is.list(U)){
      # Only subsets of variables given by the user in U, transform the list into compatible matrix
      U <- t(sapply(U, function(u){
        row <- rep(0, p)
        row[u] <- 1
        return(row)
        }))
    }else{
    # We have to compute Sobol first-order or total index
      order <- U
      if (order==0){
        # Total index
        U <- matrix(1,p,p)
        diag(U) <- 0
      }else{
        if (order==1){
          # First-order index
          U <- diag(rep(1,p))
        }else{
          stop("U is neither null, a matrix, a list, or equal to 0 or 1")
        }
      }
    }
  }
  Ul <- matrix(as.logical(U),nrow=nrow(U))
  
  # Compute sensivitiy indices for all subsets of variables defined in U
  # --------------------------------------------------------------------
  S <- apply(Ul,1,function(id){
    # id corresponds to the vector of input variables (subset) which appear in the Sobol index we compute
    # id.var contains their corresponding column indices in the X matrix 
    # (useful for categorical variables which are expanded on multiple columns via one-hot encoding)
    id.cols <- unlist(id.var[id])
    if (method=="rank"){
      # Use ranking or TSP
      if (length(id.cols)>1){
        # Use TSP
        d <- as.matrix(dist(X[,id.cols,drop=FALSE]))
        tsp <- TSP::TSP(d)
        tour <- as.integer(TSP::solve_TSP(tsp))
        id.sort <- cbind(tour,tour[c(2:n,1)])
      }else{
        # Use ranking
        id.sort <- sort(X[,id.cols], index.return = TRUE)$ix
        id.sort <- matrix(c(id.sort,id.sort[c(2:n,1)]),ncol=2)
      }
      n.knn <- 2
    }else{
      # Use knn
      if (n<=n.limit){
        # Direct search of nearest neighbors
        d <- as.matrix(dist(X[,id.cols,drop=FALSE]))
        id.sort <- t(apply(d,2,order)[1:n.knn,])
      }else{
        # Approximate search of nearest neighbors
        id.sort <- RANN::nn2(X[,id.cols,drop=FALSE],k=n.knn)$nn.idx
      }
    }
    return(1-meanknn(Y,id.sort,n.knn)) # Estimate of Var(E(Y|Xu))/Var(Y)
  }) 

  # Compute Shapley effects if asked, or Sobol first-order/total
  # --------------------------------------------------------------
  if (return.shap){
    if (!randperm){
      # Compute Shapley effects for each input variable by using all the subsets Sobol indices computed above
      Su <- matrix(0,q,2^p)
      Su[,2:2^p] <- S
      if (!noise){
        # Replace estimate of Var(E(Y|X1,...,Xp))/Var(Y) by 1
        Su[,2^p] <- 1
      }
      Uu <- rbind(rep(0,p),U)
      Shap <- matrix(NA,q,p)
      for (i in 1:p){
        ids <- matrix(which(Uu[,i]==0)) # i is not in u
        diffSu <- apply(ids,1,function(id){
          U.id <- Uu[id,]
          U.plus <- U.id
          U.plus[i] <- 1 # u with i added
          id.plus <- which(apply(Uu,1,function(x) all(x==U.plus)))
          return((Su[,id.plus]-Su[,id])/choose(p-1,sum(U.id)))
        })
        Shap[,i]=apply(matrix(diffSu,nrow=q),1,sum)/p
      }
    }
    else{
      # Compute Shapley effects for each input variable by using the random permutation subsets
      Su <- cbind(rep(0,q),matrix(S,nrow=q))
      Uu <- rbind(rep(0,p),U)
      if (!noise){
        # Replace estimate of Var(E(Y|X1,...,Xp))/Var(Y) by 1
        id.noise <- which(apply(Uu,1,function(x) all(x==matrix(1,1,p))))
        Su[,id.noise] <- 1
      }
      Shap <- matrix(NA,q,p)
      for (i in 1:p){
        diffSu <- apply(perms,1,function(perm){
          idi <- which(perm==i)
          if (idi==1){
            U.id <- rep(0,p)
          }else{
            U.id <- rep(0,p)
            U.id[perm[1:(idi-1)]] <- 1
          }
          id.pi.perm <- which(apply(Uu,1,function(x) all(x==U.id)))
          U.plus <- rep(0,p)
          U.plus[perm[1:idi]] <- 1
          id.plus <- which(apply(Uu,1,function(x) all(x==U.plus)))
          return(Su[,id.plus]-Su[,id.pi.perm])
        })
        Shap[,i]=apply(matrix(diffSu,nrow=q),1,sum)/n.perm
      }
    }
    return(list(U=U,S=S,Shap=Shap,order=order))
  }else{
    # No Shapley effet computation, just return Sobol indices
    if (is.null(order)){
      # All Sobol indices for subsets of variables given by the user in U
      return(list(U=U,S=S,order=order))
    }else{
      if (order==0){
        # Total index
        S <- 1 - S
        return(list(S=S,order=order))
      }else{
        # First-order index
        return(list(S=S,order=order))
      }
    }
  }
}

tell.sobolshap_knn <- function(x, y = NULL, ...) {
  
  id <- deparse(substitute(x))
  if (! is.null(y)) {
    x$y <- y
  } 
  else if (is.null(x$y)) {
    stop("y not found")
  }
  
  n <- nrow(x$X)
  p <- ncol(x$X)
  if (is.null(dim(x$y))){
    x$y <- matrix(x$y,ncol=1)
  }
  q <- ncol(x$y)
  
  if (is.null(colnames(x$y))){
    colnames(x$y) <- paste0("Y",1:q)
  }
  
  data <-cbind(x$X,x$y)
  res <- estim.sobolshap_knn(data, 1:n, q=q, id.cat=x$id.cat, U=x$U, method=x$method, n.knn=x$n.knn,
                             return.shap=x$return.shap, randperm=x$randperm, n.perm=x$n.perm, rescale=x$rescale, 
                             n.limit=x$n.limit, noise=x$noise)
  
  x$S <- matrix(res$S,nrow=q)
  if (!is.null(res$order)){
    colnames(x$S) <- colnames(x$X)
    if (q > 1){
      rownames(x$S) <- colnames(x$y)
    }
  }else{
    colnames(x$S) <- paste0("U",1:nrow(res$U))
    if (q > 1){
      rownames(x$S) <- colnames(x$y)
    }
    x$U <- res$U
    rownames(x$U) <- paste0("U",1:nrow(res$U))
    colnames(x$U) <- colnames(x$X)
    if (x$return.shap){
      x$Shap <- matrix(res$Shap,nrow=q)
      colnames(x$Shap) <- colnames(x$X)
      if (q > 1){
        rownames(x$Shap) <- colnames(x$y)
      }
    }
  }
  x$order <- res$order
  
  assign(id, x, parent.frame())
  return(x)
}

extract.sobolshap_knn <- function(x, ...){
  
  if (!x$return.shap & !x$randperm){
    stop("Can only extract first-order and total indices if all combinations have been computed with return.shap=TRUE")
  }
  
  U <- x$U
  S <- x$S
  p <- ncol(U)
  q <- nrow(S)
  S1 <- matrix(NA,q,p)
  ST <- matrix(NA,q,p)
  U.idall <- which(apply(U,1,function(x) all(x==rep(1,p)))) # if noise was FALSE, the associated index will be equal to 1
  for (i in 1:p){
    # First-order
    # Get combination where only i appears 
    target.ids <- matrix(0,1,p)
    target.ids[i] <- 1
    id <- which(apply(U,1,function(x) all(x==target.ids)))
    S1[,i] <- S[,id]
    # Total effect
    # Get combination where i does not appear but all others do
    target.ids <- matrix(1,1,p)
    target.ids[i] <- 0
    id <- which(apply(U,1,function(x) all(x==target.ids)))
    ST[,i] <- S[,U.idall] - S[,id]
  }
  colnames(S1) <- colnames(ST) <- colnames(U)
  rownames(S1) <- rownames(ST) <- rownames(S)
  return(list(S=S1,T=ST))
}

print.sobolshap_knn <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (! is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    if (! is.null(x$order)) {
      cat("\n\n\n",switch(x$order+1,"Total indices:","First-order indices:"),"\n")
      print(x$S)
    }else{
      if (!is.null(x$Shap)){
        cat("\n\n\nShapley effects: \n")
        print(x$Shap)
      }else{
        cat("(empty)\n")
      }
    }
  }else{
    cat("(empty)\n")
  }
}

plot.sobolshap_knn <- function(x, ylim = c(0, 1), type.multout="lines", ...) {
  
  if (! is.null(x$y)) {
    if (!is.null(x$order)){
      q <- nrow(x$S)
      if (q==1){
        s <- as.data.frame(t(x$S)); colnames(s) <- "original"
        nodeplot(s, ylim = ylim)
        legend(x = "topright", legend = switch(x$order+1,"Total indices","First-order indices"))
      }else{
        p <- ncol(x$S)
        plot(0,ylim=c(0,1),xlim=c(1,q),main=switch(x$order+1,"Total indices","First-order indices"),ylab="",xlab="",type="n",xaxt = "n")
        axis(1, at=1:q, rownames(x$S))
        for (i in 1:p){
          if (type.multout=="points"){
            points(x$S[,i],col=i,xaxt = "n")
          }
          if (type.multout=="lines"){
            lines(x$S[,i],col=i,xaxt = "n")
          }
        }
        legend(x = "topright", legend = colnames(x$S), lty=1, col=1:p, cex=0.6)
      }
    }else{
      if (!is.null(x$Shap)){
        q <- nrow(x$Shap)
        if (q==1){
          shap <- as.data.frame(t(x$Shap)); colnames(shap) <- "original"
          nodeplot(shap, ylim = ylim)
          legend(x = "topright", legend = "Shapley effects")
        }else{
          p <- ncol(x$Shap)
          plot(0,ylim=c(0,1),xlim=c(1,q),main="Shapley effects",ylab="",xlab="",type="n",xaxt = "n")
          axis(1, at=1:q, rownames(x$Shap))
          for (i in 1:p){
            if (type.multout=="points"){
              points(x$Shap[,i],col=i,xaxt = "n")
            }else{
              lines(x$Shap[,i],col=i,xaxt = "n")
            }
          }
          legend(x = "topright", legend = colnames(x$Shap), lty=1, col=1:p, cex=0.6)
        }
      }else{
        stop("No plot available")
      }
    }
  }
}

ggplot.sobolshap_knn <- function(x, ylim = c(0, 1), type.multout="lines", ...) {
  
  if (! is.null(x$y)) {
    if (!is.null(x$order)){
      q <- nrow(x$S)
      if (q==1){
        s <- as.data.frame(t(x$S)); colnames(s) <- "original"
        nodeggplot(list(s), xname = switch(x$order+1,"Total indices","First-order indices"), ylim = ylim)
      }else{
        if (q<=10){
          angle <- 0
          hjust <- 0
        }
        if (q>10 & q <=20){
          angle <- 45
          hjust <- 1
        }
        if (q>20){
          angle <- 90
          hjust <- 1
        }
        p <- ncol(x$S)
        #        df <- cbind(reshape2::melt(as.data.frame(x$S),measure.vars = 1:p),output=rep(1:q,p))
        output <- rep(1:q,p)
        df <- cbind(reshape2::melt(as.data.frame(x$S),measure.vars = 1:p),output)
        if (type.multout=="points"){
          ggplot(df,aes(x=output,y=value,color=variable)) + geom_point() + 
            coord_cartesian(ylim=ylim) + labs(y="", x = "", title = switch(x$order+1,"Total indices","First-order indices")) + theme_bw() + 
            scale_x_continuous(labels = rownames(x$S), breaks = 1:q) +
            theme(axis.text.x = element_text(angle = angle, hjust = hjust))
        }else{
          ggplot(df,aes(x=output,y=value,color=variable)) + geom_line() +
            coord_cartesian(ylim=ylim) + labs(y="", x = "", title = switch(x$order+1,"Total indices","First-order indices")) + theme_bw() + 
            scale_x_continuous(labels = rownames(x$S), breaks = 1:q) +
            theme(axis.text.x = element_text(angle = angle, hjust = hjust))
          
        }
      }
    }else{
      if (!is.null(x$Shap)){
        q <- nrow(x$Shap)
        if (q==1){
          shap <- as.data.frame(t(x$Shap)); colnames(shap) <- "original"
          nodeggplot(list(shap), xname = "Shapley effects", ylim = ylim)
        }else{
          if (q<=10){
            angle <- 0
            hjust <- 0
          }
          if (q>10 & q <=20){
            angle <- 45
            hjust <- 1
          }
          if (q>20){
            angle <- 90
            hjust <- 1
          }
          p <- ncol(x$Shap)
          output <- rep(1:q,p)
          df <- cbind(reshape2::melt(as.data.frame(x$Shap),measure.vars = 1:p),output)
          if (type.multout=="points"){
            ggplot(df,aes(x=output,y=value,color=variable)) + geom_point() + 
              coord_cartesian(ylim=ylim) + labs(y="", x = "", title = "Shapley effects") + theme_bw() + 
              scale_x_continuous(labels = rownames(x$Shap), breaks = 1:q) +
              theme(axis.text.x = element_text(angle = angle, hjust = hjust))
          }else{
            ggplot(df,aes(x=output,y=value,color=variable)) + geom_line() +
              coord_cartesian(ylim=ylim) + labs(y="", x = "", title = "Shapley effects") + theme_bw()+ 
              scale_x_continuous(labels = rownames(x$Shap), breaks = 1:q) +
              theme(axis.text.x = element_text(angle = angle, hjust = hjust))
          }
        }
      }else{
        stop("No plot available")
      }
    }
  }
}