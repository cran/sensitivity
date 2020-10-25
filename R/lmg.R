# Lindeman-Merenda-Gold indices
#
# Bertrand Iooss 2020

# somme sur permutations avec r2
estim1.lmg <- function(data, logistic, i = 1:nrow(data) ) { # logistic sert a rien  
  d <- data[i, ]
  p <- ncol(d) - 1
  lmg <- numeric(p)
  
  perms <- gtools::permutations(p, p, 1:p) # Generate all d! permutations
  m <- nrow(perms)
  
  for (j in 1:p) {
    
    for (i in 1:m) {
      
      ind <- p + 1
      k <- 1
      while (perms[i,k] != j){
        ind <- c(ind,k)
        k <- k + 1
      }
      Xtildej.lab <- paste(colnames(d[c(-1, -ind, -j-1)]), collapse = "+")
      #      if (Xtildej.lab == "") Xtildej.lab <- "."
      
      if (Xtildej.lab != ""){
        lm.Xj <- lm(formula(paste(colnames(d)[j+1], "~", Xtildej.lab)), data = d)
        lmg[j] <- lmg[j] + cor(d[1], d[j+1] - fitted(lm.Xj))^2
      }
    }
  }
  lmg / factorial(p)
}

# somme sur combinaisons avec r2
estim.lmg <- function(data, logistic, i = 1:nrow(data) ) { # logistic sert a rien  
  d <- data[i, ]
  p <- ncol(d) - 1
  lmg <- numeric(p)
  
  for (j in 1:p) {
    
    for (i in 0:(p-1)) {
     
      comb <- combn(c(1:p)[-j], i) # Generate all combinations without j
      m <- ncol(comb)
      
      for (k in 1:m) {
        Xtildej.lab <- paste(colnames(d[comb[,k] + 1]), collapse = "+")
        if (Xtildej.lab == "") Xtildej.lab <- "."
        
        lm.Xj <- lm(formula(paste(colnames(d)[j+1], "~", Xtildej.lab)), data = d)
      
        lmg[j] <- lmg[j] + factorial(p-1-i) * factorial(i) * cor(d[1], d[j+1] - fitted(lm.Xj))^2
      }
    }
  }
  lmg / factorial(p)
}


lmg <- function(X, y, logistic = FALSE, nboot = 0, conf = 0.95) {
  data <- cbind(Y = y, X)
  
  if (nboot == 0) {
    lmg <- data.frame(original = estim.lmg(data, logistic))
    rownames(lmg) <- colnames(X)
  } else {
    boot.lmg <- boot(data, estim.lmg, logistic = logistic, R = nboot)
    pcc <- bootstats(boot.lmg, conf, "basic")
    rownames(lmg) <- colnames(X)
  }

  out <- list(X = X, y = y, nboot = nboot, conf = conf,
              call = match.call())
  class(out) <- "lmg"
  out$LMG <- lmg
  
  return(out)
}


print.lmg <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if ("LMG" %in% names(x)) {
    cat("\nLindeman-Merenda-Gold indices (LMG):\n")
    print(x$LMG)
  }
}


plot.lmg <- function(x, ylim = c(-1,1), ...) {  
  if ("LMG" %in% names(x)) {
    nodeplot(x$LMG, ylim = ylim, main = "LMG")
  }
}


ggplot.lmg <- function(x, ylim = c(-1,1), ...) {  
  if ("LMG" %in% names(x)) {
    nodeggplot(listx = list(x$PCC), xname = "LMG", ylim = ylim, title = "LMG")
  }
}

