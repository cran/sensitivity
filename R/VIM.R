##############################################
# Variance-based Importance Measures in linear regression model 
# for data analysis context
#
# Bertrand Iooss 2026
############################################
# Needed packages
#library(car)
#library(boot)
#library(parallel)


VIM <- function(X, y, logistic = FALSE, nboot = 0, conf = 0.95, max.iter=1000, parl=NULL){
  
  dat <- data.frame(y,X)
  
  if (!logistic){
    modlin <- glm(y~., data=dat)
  } else{
    fit.control=glm.control(maxit=max.iter)
    modlin <- glm(y~., data = dat, family="binomial", control=fit.control)
  }
   
  R2 <- 1 - (modlin$deviance/modlin$null.deviance)
  Q2 = 1 - boot::cv.glm(dat,modlin)$delta[1]/var(y)
  
  vif <- car::vif(modlin)
  src2 <- src(X = X, y = y, logistic = logistic, nboot = nboot, conf = conf)$SRC^2
  pcc2 <- pcc(X = X, y = y, logistic = logistic, nboot = nboot, conf = conf)$PCC^2
  lmg <- lmg(X = X, y = y, logistic = logistic, nboot = nboot, conf = conf, 
             max.iter = max.iter, parl = parl)
  pmvd <- pmvd(X = X, y = y, logistic = logistic, nboot = nboot, conf = conf, 
               max.iter = max.iter, parl = parl)
  

  out <- list(X = X, y = y, logistic = logistic, nboot = nboot, conf = conf, 
              max.iter = max.iter, parl = parl, 
              call = match.call())
  class(out) <- "VIM"
  out$R2 <- R2
  out$Q2 <- Q2
  out$VIF <- vif
  out$SRC2 <- src2
  out$PCC2 <- pcc2
  if (nboot == 0){
    out$LMG <- lmg$lmg
    out$PMVD <- pmvd$pmvd
  } else{
    out$LMG <- lmg$conf_int
    out$PMVD <- pmvd$conf_int
  }
  
  return(out)
}


print.VIM <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (!is.null(x$R2)) cat("\nDetermination coefficient of linear regression model: R2 = ", x$R2, "\n")
  if (!is.null(x$Q2)) cat("\nPredictivity coefficient of linear regression model: Q2 = ", x$Q2, "\n")
  if ("VIF" %in% names(x)) {
    cat("\nVariance Inflation Factors (VIF):\n")
    print(x$VIF)
  }
  if ("SRC2" %in% names(x)) {
    cat("\nSquared Standard Regression Coefficients (SRC2):\n")
    print(x$SRC2)
  }
  if ("PCC2" %in% names(x)) {
    cat("\nSquared Partial Correlation Coefficients (PCC2):\n")
    print(x$PCC2)
  }
  if ("LMG" %in% names(x)) {
    cat("\nLindeman, Merenda and Gold indices (LMG):\n")
    print(x$LMG)
  }
  if ("PMVD" %in% names(x)) {
    cat("\nProportional Marginal Variance Decomposition indices (PMVD):\n")
    print(x$PMVD)
  }
}


plot.VIM <- function(x, ylim = c(0, 1), ...) {
  p <- ncol(x$X)
  pch = c(21, 24)
  nodeplot(as.data.frame(x$LMG), xlim = c(1, p + 1), ylim = ylim, pch = pch[1])
  nodeplot(as.data.frame(x$PMVD), xlim = c(1, p + 1), ylim = ylim, labels = FALSE,
             pch = pch[2], at = (1:p)+.3, add = TRUE)
  legend(x = "topright", legend = c("LMG", "PMVD"), pch = pch)
}

ggplot.VIM <- function(data, mapping = aes(), ..., ylim = c(0,1), environment = parent.frame()) {
  x <- data
  p <- ncol(x$X)
  pch = c(21, 24)
  nodeggplot(listx = list(as.data.frame(x$LMG),as.data.frame(x$PMVD)), xname = c("LMG","PMVD"), ylim = ylim, pch = pch)
}
