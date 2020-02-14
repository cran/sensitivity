# compute squaredCoef i generalized chaos
# = < y, e_{1, l1}...e_{d, ld}>^2 
# = [1/lambda_{i,li} < dy/dxi, e_{1,l1}... e'_{i,li} ... e_{d,ld}>]^2
#
# Authors: Olivier Roustant and Bertrand Iooss (2019)
#
# Reference: 
# O. Roustant, F. Gamboa and B. Iooss, Parseval inequalities and lower bounds 
# for variance-based sensitivity indices, Preprint, ARXIV: 1906.09883

PoincareChaosSqCoef <- function(PoincareEigen, multiIndex, 
                                design, output, outputGrad = NULL, 
                                inputIndex = 1, der = FALSE,
                                method = "unbiased"){
  # PoincareEigen: output list from PoincareOptimal
  # multiIndex: vector of indices (l1, ..., ld)
  # design: design of experiments (matrix of size n x d)
  # output: vector of length n (y1, ..., yn) of observations at design = (x1, ..., xn)
  # outputGrad: matrix n x d whose columns contain the partial derivatives at X
  # inputIndex: index of the input variable (between 1 and d)
  # der: should we use the formula with derivatives to compute the square coefficient ?
  d <- length(PoincareEigen)
  multiIndex <- as.integer(multiIndex)
  if (length(multiIndex) != d) stop("The length of multiindex must be equal to the number of input variables")
  nonZeroSet <- which(multiIndex != 0)
  
  if (multiIndex[inputIndex] == 0 & der) stop("Division by zero. Change multiIndex[inputIndex] or der values")
  
  chaos <- 1
  if (!der){
    for (i in nonZeroSet){
      basis <- approxfun(x = PoincareEigen[[i]]$knots, 
                         y = PoincareEigen[[i]]$vectors[, multiIndex[i] + 1])
      chaos <- chaos * basis(design[, i])
    }
    res <-  output * chaos
  } else {
    if (multiIndex[inputIndex] == 0) {   
      res <- 0   # the derivative of the constant (one) eigenvalue is zero
    } else {
      for (i in setdiff(nonZeroSet, inputIndex)){
        basis <- approxfun(x = PoincareEigen[[i]]$knots, 
                           y = PoincareEigen[[i]]$vectors[, multiIndex[i] + 1])
        chaos <- chaos * basis(design[, i])
      }
      i <- inputIndex
      basis <- approxfun(x = PoincareEigen[[i]]$knots, 
                         y = PoincareEigen[[i]]$der[, multiIndex[i] + 1])
      chaos <- chaos * basis(design[, i]) / PoincareEigen[[i]]$values[multiIndex[i] + 1]
      res <- outputGrad[, i] * chaos
    }
  }
  
  c2 <- squaredIntEstim(res, method = method)
  return(c2)
}
