linsa <- function(model = NULL, x, pcc = TRUE, rank = FALSE,
                  nboot = 0, conf = 0.95, ...)
{
  # OBJECT OF CLASS "linsa"
  
  sa <- list(model = model, x = x, rank = rank, nboot = nboot, conf = conf,
             y = NULL, src = NULL, call = match.call())
  if (pcc == TRUE){
    sa <- c(sa, list(pcc = NULL))
  }
  class(sa) <- "linsa"
  
  # COMPUTATION OF THE RESPONSE AND SENSITIVITY ANALYSIS
  
  if (!is.null(sa$model)){
      response(sa, ...)
      compute(sa)
  }
  
  # RETURN OF THE OBJECT OF CLASS "linsa"
  
  return(sa)
}


## estim.pear <- function(data, i)
## {
##   d <- data[i, ]
##   p <- ncol(d) - 1
##   out <- numeric(p)
##   for (j in 1 : p)
##     out[j] <- cor(d[1], d[j+1])
##   return(out)
## }


estim.src <- function(data, i)
{
  d <- data[i, ]
  f <- formula(paste(colnames(d)[1], "~", paste(colnames(d)[-1],
                                                collapse = "+")))
  reg <- lm(formula = f, data = d)
  return(coefficients(reg)[-1] * sd(d[-1]) / sd(d[1]))
}


estim.pcc <- function(data, i)
{
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


compute.linsa <- function(sa, y = NULL)
{
  id <- deparse(substitute(sa))

  # EXTERNAL MODEL

  if (! is.null(y))
    sa$y <- y
  
  # RANK ANALYSIS ?

  if (sa$rank == TRUE){
    for (i in 1:ncol(sa$x))
      sa$x[,i] <- rank(sa$x[,i])
    sa$y <- rank(sa$y)
  }

  # ESTIMATION OF THE INDICES
  
  data <- as.data.frame(cbind(Y=sa$y, sa$x))
##   sa$pear <- estim(estim.pear, data, sa$nboot, sa$conf)
##   rownames(sa$pear) <- colnames(x)
  sa$src <- estim(estim.src, data, sa$nboot, sa$conf)
  rownames(sa$src) <- colnames(sa$x)
  if ("pcc" %in% names(sa)){
    sa$pcc <- estim(estim.pcc, data, sa$nboot, sa$conf)
    rownames(sa$pcc) <- colnames(sa$x)
  }
  
  # SAVING OF THE INDICES
  
  assign(id, sa, parent.frame())
}


print.linsa <- function(x, ...)
{
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


plot.linsa <- function(x, ask = TRUE, ...)
{
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
