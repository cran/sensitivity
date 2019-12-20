squaredIntEstim <- function(x, method = "unbiased"){
  res <- (mean(x))^2
  if (method == "unbiased"){  
    n <- length(x)
    res <- res - var(x)/n
  }
  return(res)
}
