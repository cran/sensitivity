# GENERIC METHOD FOR SENSITIVITY ANALYSIS

compute <- function(sa, y = NULL)
  UseMethod("compute")


# COMPUTING THE RESPONSE OF A 'SA' OBJECT

response <- function(sa, ...)
{
  id <- deparse(substitute(sa))
  
  # old version :
##   if (class(sa$model) == "function")
##     y <- sa$model(sa$x, ...)
##   else
##     y <- predict(sa$model, sa$x, ...)

  # new version :
  if (class(sa$model) == "function")
     y <- sa$model(sa$x, ...)
   else if (TRUE %in% (paste("predict.", class(sa$model), sep="") %in%
                       methods(predict)))
     y <- predict(sa$model, sa$x, ...)
   else
     stop("The model isn't a function or a predictor")

  if (class(y) != "numeric"){
    y <- as.numeric(y)
    warning("Conversion of the response to numeric")
  }

  sa$y <- y
  assign(id, sa, parent.frame())
}


# BOOTSTRAP STATISTICS
#
# b : object of class 'boostrap'
# confidence : confidence level for bootstrap bias-corrected confidence
#   intervals
# type : type of confidence interval, "norm" or "basic"
#
# returns : a data.frame of bootstrap statistics

bootstats <- function(b, conf = 0.95, type = "norm")
{
  p <- length(b$t0)
  lab <- c("original", "bias", "std. error", "min. c.i.", "max. c.i.")
  out <-  as.data.frame(matrix(nrow = p, ncol = length(lab),
                               dimnames = list(NULL, lab)))

  for (i in 1:p){
    
      # original estimation, bias, standard deviation
      
      out[i, "original"] <- b$t0[i]
      out[i, "bias"] <- mean(b$t[, i]) - b$t0[i]
      out[i, "std. error"] <- sd(b$t[, i])
      
      # confidence interval

      if (type == "norm"){
        ci <- boot.ci(b, index = i, type = "norm", conf = conf)
        if (!is.null(ci)){
          out[i, "min. c.i."] <- ci$norm[2]
          out[i, "max. c.i."] <- ci$norm[3]
        }
      }
      else if (type == "basic"){
        ci <- boot.ci(b, index = i, type = "basic", conf = conf)
        if (!is.null(ci)){
          out[i, "min. c.i."] <- ci$basic[4]
          out[i, "max. c.i."] <- ci$basic[5]
        }
      }
    }
  
  return(out)
}


# ESTIMATING FUNCTION
#
# estim.fun : function(data, i, ...) with i the vector of lines for bootstraping
# data : data.frame or matrix
# nboot : nb of bootstrap replicates (0 if no bootstrap)
# conf : confidence level for confidence intervals
# ... : more arguments passed to estim.fun 
#
# returns : a data.frame of estimated statistics

estim <- function(estim.fun, data, nboot, conf=0.95, ...)
{
  if (nboot == 0){
    
    # estimation without bootstrap
    
    out <- data.frame(original = estim.fun(data, 1:nrow(data), ...))
  }
  else{

    # estimation with bootstrap
    
    b <- boot(data, estim.fun, R = nboot, ...)
    out <- bootstats(b, conf, "basic")
  }

  return(out)
}


# PLOTTING CONFIDENCE INTERVALS

nodeplot <- function(x, xlim = NULL, ylim = NULL, labels = TRUE,
                     col = par("col"), pch = 21, bg = "white",
                     add = FALSE, at = NULL)
{
  n <- nrow(x)
  if (is.null(xlim))
    xlim <- c(1, n)
  if (is.null(ylim))
    ylim <- c(min(x), max(x))
  if (is.null(at))
    at <- 1 : n
  if (add)
    par(new = TRUE)

  # axes
  
  plot(0, xlim = xlim, ylim = ylim, axes = FALSE,
       xlab = "", ylab = "", main = "", type = "n")
  if (class(labels) == "logical")
    if (labels == TRUE)
      axis(side = 1, at = at, labels = rownames(x))
    else
      axis(side = 1, at = at, labels = FALSE, tick = FALSE)
  else if (class(labels) == "character")
    axis(side = 1, at = at, labels = labels)
  axis(side = 2)
  box()

  # to bias or not to bias...

  if ("bias" %in% colnames(x))
    xx <- x[["original"]] - x[["bias"]]
  else
    xx <- x[["original"]]

  # confidence intervals

  if (("min. c.i." %in% colnames(x)) & "max. c.i." %in% colnames(x))
    for (i in 1 : n)
      lines(c(at[i], at[i]), c(x[["min. c.i."]][i], x[["max. c.i"]][i]),
            col = col)

  # points

  points(at, xx, col = col, pch = pch, bg = bg)
}


# PLOTTING CONFIDENCE INTERVALS FOR 2 STATISTICS

crossplot <- function(x, y, xlim = NULL, ylim = NULL, xlab = "", ylab = "",
                     col = par("col"), pch = 21, bg = "white",
                     add = FALSE, labels = NULL)
{
  n <- nrow(x)
  if (is.null(xlim))
    xlim <- c(min(x), max(x))
  if (is.null(ylim))
    ylim <- c(min(y), max(y))
  if (add)
    par(new = TRUE)

  # axes
  
  plot(0, xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab, main = "", type = "n")

  # to bias or not to bias...

  if ("bias" %in% colnames(x)){
    xx <- x[["original"]] - x[["bias"]]
    yy <- y[["original"]] - y[["bias"]]
  }
  else{
    xx <- x[["original"]]
    yy <- y[["original"]]
  }

  # confidence intervals
  
  if (("min. c.i." %in% colnames(x)) & "max. c.i." %in% colnames(x))
    for (i in 1 : n)
      lines(c(x[["min. c.i."]][i], x[["max. c.i"]][i]), c(yy[i], yy[i]),
            col = col)
  if (("min. c.i." %in% colnames(y)) & "max. c.i." %in% colnames(y))
    for (i in 1 : n)
      lines(c(xx[i], xx[i]), c(y[["min. c.i."]][i], y[["max. c.i"]][i]),
            col = col)

  # points
  
  if (! is.null(labels))
    text(xx, yy, labels = labels)
  else
    points(xx, yy, col = col, pch = pch, bg = bg)
}
