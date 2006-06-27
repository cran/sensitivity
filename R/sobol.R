subsets <- function(set, size, type = "eq") # finds all the subsets of
{                                           # size (less or) equal to n
  n <- length(set)
  if ( size < 1 | size > n )
    stop()
  else if ( size == 1 )
    return(as.list(set))
  else if ( type == "leq" )
    return(c(Recall(set, size - 1, type = "leq"), Recall(set, size)))  
  else if ( size == n )
    return(list(set))
  else{
    out <- list()
    for ( i in 1 : (n - size + 1) )
      for ( j in Recall(set[(i + 1) : n], size - 1) )
        out <- c(out, list(c(set[i], j)))
    return(out)
  }
}


sobol <- function(model = NULL, x1, x2, order = 1, nboot = 0, conf = 0.95, ...)
{
  # ARGUMENTS

  if ( (ncol(x1) != ncol(x2)) | (nrow(x1) != nrow(x2)) )
    error("The two samples x1 and x2 must have the same dimensions")
  p <- ncol(x1)

  # DESIGN OF EXPERIMENTS
  
  x <- x1
  for ( i in subsets(set = 1 : p, size = order, type = "leq") ){
    xb <- x2
    xb[, i] <- x1[, i]
    x <- rbind(x, xb) 
  }

  # OBJECT OF CLASS "sobol"
  
  sa <- list(model = model, x1 = x1, x2 = x2, order = order, nboot = nboot,
             conf = conf, , x = x, y = NULL, S = NULL, call = match.call())

  class(sa) <- "sobol"

  # COMPUTATION OF THE RESPONSE AND SENSITIVITY ANALYSIS
  
  if ( !is.null(sa$model) ){
      response(sa, ...)
      compute(sa)
  }

  # RETURN OF THE OBJECT OF CLASS "sobol"
  
  return(sa)
}


estim.sobol <- function(data, i, mask)
{
  # 1. computation of subset indices
  d <- as.matrix(data[i, ]) # as.matrix pour que colSums renvoie un numeric...
  n <- nrow(d)
  V <- var(d[, 1])
  m2 <- mean(d[, 1])^2
  subset.indices <- (colSums(d[, -1] * d[, 1]) / (n - 1) - m2) / V

  # 2. computation of multi-order indices, based on subset ones
  multi.indices <- subset.indices
  ni <- length(multi.indices)
  for ( j in 1 : ni )
    multi.indices[j] <- sum(multi.indices * mask[[j]])
  
  return(multi.indices)
}


compute.sobol <- function(sa, y = NULL)
{
  id <- deparse(substitute(sa))
  
  # EXTERNAL MODEL

  if ( ! is.null(y) )
    sa$y <- y

  # ESTIMATION OF THE INDICES
  
  # indices to compute
  list.indices <- subsets(set = 1 : ncol(sa$x1), size = sa$order, type = "leq")
  ni <- length(list.indices)

  # mask used to compute hight-order indices from the indices of lower order
  mask <- list()
  for ( i in list.indices ){
    v <- rep(0, ni)
    v[match(list(i), list.indices)] <- 1
    m <- length(i)
    if ( m > 1 )
      for ( j in subsets(set = i, size = m - 1, type = "leq") )
        v[match(list(j), list.indices)] <- -1
    mask <- c(mask, list(v))
  }

  # labels
  labels.indices <- character(ni)
  factor.names <- colnames(sa$x1)
  for ( i in 1 : ni )
    labels.indices[i] <- paste(factor.names[list.indices[[i]]], collapse = ",")
 
  # estimation
  data <- matrix(sa$y, nr = nrow(sa$x1), nc = ni + 1)
  sa$S <- estim(estim.sobol, data, sa$nboot, sa$conf, mask = mask)
  rownames(sa$S) <- labels.indices
    
  # SAVING OF THE INDICES
  
  assign(id, sa, parent.frame())
}


print.sobol <- function(x, ...)
{
  cat("\nSOBOL METHOD\n")
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  cat("\nModel runs:", length(x$y), "\n")
  cat("\nAnova decomposition indices:\n")
  print(x$S)
}


plot.sobol <- function(x, ...)
{
  nodeplot(x$S, ylim = c(0, 1))
}
