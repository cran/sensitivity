# Base functions of the sensitivity package
# Gilles Pujol, 2006 - 2008


ask <- function(x, ...)
  UseMethod("ask")


tell <- function(x, y = NULL, ...)
  UseMethod("tell")


response <- function(x, loop = FALSE, ...) {
  id <- deparse(substitute(x))
  
  if (class(x$model) == "function") {
    if (loop) {
      n <- nrow(x$X)
      y <- numeric(n)
      for (i in 1:n) {
        y[i] <- x$model(x$X[i,], ...)
      }
    } else {
      y <- x$model(x$X, ...)
    }
  } else if (TRUE %in% (paste("predict.", class(x$model), sep="") %in% methods(predict))) {
    y <- predict(x$model, x$X, ...)
  } else {
    stop("The model isn't a function or does not have a predict method")
  }

  if (class(y) != "numeric") {
    y <- as.numeric(y)
    warning("Conversion of the response to numeric")
  }

  x$y <- y
  assign(id, x, parent.frame())
}
