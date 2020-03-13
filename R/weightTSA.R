# Author : Bertrand Iooss (2020)

weightTSA <- function(Y, c, upper = TRUE, type = "indicTh", param = 1) {
  # Y = the output vector
  # c = the threshold
  # upper = TRUE for upper threshold and FALSE for lower threshold
  # type = the weight function type ("indicTh", "zeroTh", logistic", "exp1side")
  # indicTh: indicator-thresholding, zeroTh: zero-thresholding,
  # logistic: logistic transf. at the threshold,
  # exp1side: exponential transf. below the threshold
  # param = the parameter value for "logistic" and "exp1side" types

  if (upper){
    if (type == "indicTh") wY <- as.numeric(Y>c)
    if (type == "zeroTh") wY <- Y * (Y>c)
    if (type == "logistic") wY <- 1 / (1 + exp(-param * (Y-c) / abs(c)) ) # Spagnol & Da Veiga
    if (type == "exp1side") wY <- exp( - (c-Y)*((c-Y)>0) / (param * sd(Y)/5) ) # Raguet & Marrel
  } else{
    if (type == "indicTh") wY <- as.numeric(Y<c)
    if (type == "zeroTh") wY <- Y * (Y<c)
    if (type == "logistic") wY <- 1 / (1 + exp(-param * (c-Y) / abs(c)) ) # Spagnol & Da Veiga
    if (type == "exp1side") wY <- exp( - (Y-c)*((Y-c)>0) / (param * sd(Y)/5) ) # Raguet & Marrel
  }
  return(wY)
}
