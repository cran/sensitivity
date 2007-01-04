
## File: srcpcc.R
## Description: linear sensitivity analysis
## Author: Gilles Pujol


srcpcc <- function(model = NULL, x, pcc = TRUE, rank = FALSE,
                  nboot = 0, conf = 0.95, ...){
  
  ## OBJECT OF CLASS "srcpcc"
  
  sa <- list(model = model, x = x, rank = rank, nboot = nboot, conf = conf,
             y = NULL, src = NULL, call = match.call())
  if (pcc == TRUE){
    sa <- c(sa, list(pcc = NULL))
  }
  class(sa) <- "srcpcc"
  
  ## COMPUTATION OF THE RESPONSE AND SENSITIVITY ANALYSIS
  
  if (!is.null(sa$model)){
      response(sa, ...)
      tell(sa)
  }
  
  ## RETURN OF THE OBJECT OF CLASS "srcpcc"
  
  return(sa)
}


estim.src <- function(data, i = 1 : nrow(data)){
  d <- data[i, ]
  f <- formula(paste(colnames(d)[1], "~",
                     paste(colnames(d)[-1], collapse = "+")))
  reg <- lm(formula = f, data = d)
  return(coefficients(reg)[-1] * sd(d[-1]) / sd(d[1]))
}


estim.pcc <- function(data, i = 1 : nrow(data)){
  d <- data[i, ]
  p <- ncol(d) - 1
  out <- numeric(p)
  for (j in 1 : p){
    X.mj <- paste(colnames(d)[c(-1, -j-1)], collapse = "+")
    reg1 <- lm(formula = formula(paste(colnames(d)[1], "~", X.mj)),
               data = d)
    reg2 <- lm(formula = formula(paste(colnames(d)[j+1], "~", X.mj)),
               data = d)
    out[j] <- cor(d[1] - fitted(reg1), d[j+1] - fitted(reg2))
  }
  return(out)
}


tell.srcpcc <- function(sa, y = NULL){
  id <- deparse(substitute(sa))

  ## EXTERNAL MODEL

  if (! is.null(y))
    sa$y <- y
  
  ## RANK ANALYSIS ?

  if (sa$rank == TRUE){
    for (i in 1:ncol(sa$x))
      sa$x[,i] <- rank(sa$x[,i])
    sa$y <- rank(sa$y)
  }

  ## a first linear regression (to allow the user to control the linearity)

  data <- as.data.frame(cbind(Y=sa$y, sa$x))
  
  f <- formula(paste(colnames(data)[1], "~",
                     paste(colnames(data)[-1], collapse = "+")))
  sa$lm <- lm(formula = f, data = data)

  ## ESTIMATION OF THE INDICES
  
  if (sa$nboot == 0){
    ## single estimation
    sa$src <- data.frame(original = estim.src(data))
    rownames(sa$src) <- colnames(sa$x)
    
    if ("pcc" %in% names(sa)){
      sa$pcc <- data.frame(original = estim.pcc(data))
      rownames(sa$pcc) <- colnames(sa$x)
    }
  }
  else{
    ## bootstrap estimation
    boot.src <- boot(data, estim.src, R = sa$nboot)
    sa$src <- bootstats(boot.src, sa$conf, "basic")
    rownames(sa$src) <- colnames(sa$x)
    
    if ("pcc" %in% names(sa)){
      boot.pcc <- boot(data, estim.pcc, R = sa$nboot)
      sa$pcc <- bootstats(boot.pcc, sa$conf, "basic")
      rownames(sa$pcc) <- colnames(sa$x)
    }
  }
  
  ## SAVING OF THE INDICES
  
  assign(id, sa, parent.frame())
}


print.srcpcc <- function(x, ...){
   if (x$rank == FALSE)
     titles <- c("\nStandardized Regression Coefficients:\n",
                 "\nPartial Correlation Coefficients:\n")
   else
     titles <- c("\nStandardized Rank Regression Coefficients:\n",
                 "\nPartial Rank Correlation Coefficients:\n")
  cat("\nLINEAR SENSITIVITY ANALYSIS\n")
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  cat("\nModel runs:", length(x$y), "\n")
  cat(titles[1])
  print(x$src)
  if ("pcc" %in% names(x)){
    cat(titles[2])
    print(x$pcc)
  }
}


plot.srcpcc <- function(x, ask = TRUE, ...){
  if (ask){
    op <- par(ask = TRUE)
    on.exit(par(op))
  }
  if (x$rank == FALSE)
    titles <- c("SRC", "PCC")
  else
    titles <- c("SRRC", "PRCC")
  nodeplot(x$src, ylim = c(-1, 1))
  title(main = titles[1])
  if ("pcc" %in% names(x)){
    nodeplot(x$pcc, ylim = c(-1, 1))
    title(main = titles[2])
  }
}
