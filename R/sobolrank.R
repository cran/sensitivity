sobolrank <- function(model = NULL, X, nboot = 0, conf = 0.95, nsample = round(0.8*nrow(X)), ...) {
  
  if (is.data.frame(X)){
    X <- as.matrix(unname(X))
  }else if(!is.matrix(X)){
    stop("The sample X must be a matrix or a data frame")
  }
  
  x <- list(model = model, X = X, nboot = nboot, conf = conf, nsample = nsample, call = match.call()) 
  class(x) <- "sobolrank"
  
  #calcul of the response for explicit model
  if (! is.null(x$model)) {
    response(x, ...)
    x=tell(x, ...)
  }  
  return(x)
}

estim.sobolrank <- function(data, i=1:nrow(data)) {
  n <- nrow(data)
  ptot <- ncol(data)
  p <- ptot - 1
  X <- data[i,1:p]
  Y <- data[i,ptot]
  S = matrix(0,nrow=p,ncol=1)
  EY <- mean(Y)
  VarY <- var(Y)
  
  for (j in 1:p){
    id.sort <- sort(X[,j], decreasing = FALSE, index.return = TRUE)$ix
    S[j] <- (mean(Y[id.sort]*Y[id.sort[c(2:n,1)]]) - EY^2)/VarY
  }
  return(S)
}

tell.sobolrank <- function(x, y = NULL, ...) {
  
  id <- deparse(substitute(x))
  if (! is.null(y)) {
    x$y <- y
  } 
  else if (is.null(x$y)) {
    stop("y not found")
  }
  
  n <- nrow(x$X)
  p <- ncol(x$X)
  
  data <-cbind(x$X,x$y)
  if (x$nboot == 0){
  res <- estim.sobolrank(data, 1:n)
  x$S <- data.frame(res)
  colnames(x$S) <- "original"
  rownames(x$S) <- colnames(x$X)
  }else{
    sample.boot <- function(data,mle){
      out <- data[sample.int(nrow(data),mle),,drop=FALSE]
      return(out)
    }
    S.boot <- boot(data, estim.sobolrank, R = x$nboot, sim = "parametric", ran.gen = sample.boot, mle = x$nsample)
    x$S <- bootstats(S.boot, x$conf, "basic")
  }
  
  assign(id, x, parent.frame())
  return(x)
}

print.sobolrank <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (! is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    if (! is.null(x$S)) {
      cat("\n\n\nFirst order indices: \n")
      print(x$S)
    }
    else{
      cat("(empty)\n")
    }
  }
}

plot.sobolrank <- function(x, ylim = c(0, 1), ...) {
  
  if (! is.null(x$y)) {
    nodeplot(x$S, ylim = ylim)
    legend(x = "topright", legend = "First-order indices")
  }
}

ggplot.sobolrank <- function(x, ylim = c(0, 1), ...) {
  
  if (! is.null(x$y)) {
    nodeggplot(list(x$S), xname = "First-order indices", ylim = ylim)
  }
}