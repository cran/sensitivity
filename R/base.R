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
      y <- sapply(1:n, function(i){
        x$model(x$X[i,], ...)
      })
      if(is.matrix(y)){
        # This means the model function returned vectors, so y is a 
        # (m times n)-matrix (with unknown m). For better consistency with the 
        # default case that y is a vector of length n, y is transposed to a 
        # (n times m)-matrix.
        # (Remark: Per default, sapply creates one column for every element the 
        # function is applied to.)
        y <- t(y)
      } else if(is.list(y)){
        check_matrices <- vapply(y, function(x) is.matrix(x), 
                                 FUN.VALUE = logical(1))
        if(all(check_matrices)){
          dims <- vapply(y, function(x) dim(x), FUN.VALUE = numeric(2))
          check_dims <- apply(dims, 1, function(x) length(unique(x)) == 1)
          if(all(check_dim)){
            # This means the model function returned matrices of the same sizes.
            # Let m be the number of rows and z be the number of columns.
            # For better consistency, results are transposed so that y is a list
            # of length z containing in each element a (n times m)-matrix.
            
            # Number of rows and columns of all matrices:
            dims <- dims[, 1]
            # Save the column names for restoring them later:
            y_colnames <- colnames(y[[1]])
            # y_mat is a ((n * z) times m)-matrix:
            y_mat <- t(do.call(cbind, y))
            # The rownames don't make sense anymore:
            rownames(y_mat) <- NULL
            nz <- nrow(y_mat)
            y <- lapply(1:dims[2], function(j){
              y_mat[seq(j, nz, dims[2]), , drop = FALSE]
            })
            # Restore the original column names (these are now the names of the
            # elements of y):
            names(y) <- y_colnames
          } else{
            stop("x$model returned a list of matrices, but they are not all",
                 "of the same dimensions")
          }
        } else{
          stop("x$model returned a list, but not all elements are matrices")
        }
      } else{
        stop("x$model returned an object that can't be handled")
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

########################################
# function added in 2014

plotFG <- function(x) 
  UseMethod("plotFG")

