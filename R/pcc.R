# Partial Correlation Coefficients
#
# Gilles Pujol 2006
# Bertrand Iooss 2020 for Semi-Partial Correlation Coefficients


estim.pcc <- function(data, semi, i = 1:nrow(data)) {  
  d <- data[i, ]
  p <- ncol(d) - 1
  pcc <- numeric(p)
  for (j in 1:p) {
    Xtildej.lab <- paste(colnames(d)[c(-1, -j-1)], collapse = "+")
    lm.Y <- lm(formula(paste(colnames(d)[1], "~", Xtildej.lab)), data = d)
    lm.Xj <- lm(formula(paste(colnames(d)[j+1], "~", Xtildej.lab)), data = d)
    if (! semi) {
      pcc[j] <- cor(d[1] - fitted(lm.Y), d[j+1] - fitted(lm.Xj))
    } else {
      pcc[j] <- cor(d[1], d[j+1] - fitted(lm.Xj))
    }
  }
  pcc
}


pcc <- function(X, y, rank = FALSE, semi = FALSE, nboot = 0, conf = 0.95) {
  data <- cbind(Y = y, X)

  if (rank) {
    for (i in 1:ncol(data)) {
      data[,i] <- rank(data[,i])
    }
  }
  
  if (nboot == 0) {
    pcc <- data.frame(original = estim.pcc(data, semi))
    rownames(pcc) <- colnames(X)
  } else {
    boot.pcc <- boot(data, estim.pcc, semi = semi, R = nboot)
    pcc <- bootstats(boot.pcc, conf, "basic")
    rownames(pcc) <- colnames(X)
  }

  out <- list(X = X, y = y, rank = rank, nboot = nboot, conf = conf,
              call = match.call())
  class(out) <- "pcc"
  if (! semi) {
    if (! rank) {
      out$PCC <- pcc
    } else {
      out$PRCC = pcc
    }
  } else {
    if (! rank) {
      out$SPCC <- pcc
    } else {
      out$SPRCC = pcc
    }
  }
  return(out)
}


print.pcc <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if ("PCC" %in% names(x)) {
    cat("\nPartial Correlation Coefficients (PCC):\n")
    print(x$PCC)
  } else if ("PRCC" %in% names(x)) {
    cat("\nPartial Rank Correlation Coefficients (PRCC):\n")
    print(x$PRCC)
  } else if ("SPCC" %in% names(x)) {
    cat("\nSemi-Partial Correlation Coefficients (SPCC):\n")
    print(x$SPCC)
  } else if ("SPRCC" %in% names(x)) {
    cat("\nSemi-Partial Rank Correlation Coefficients (SPRCC):\n")
    print(x$SPRCC)
  }
}


plot.pcc <- function(x, ylim = c(-1,1), ...) {  
  if ("PCC" %in% names(x)) {
    nodeplot(x$PCC, ylim = ylim, main = "PCC")
  } else if ("PRCC" %in% names(x)) {
    nodeplot(x$PRCC, ylim = ylim, main = "PRCC")
  } else if ("SPCC" %in% names(x)) {
    nodeplot(x$SPCC, ylim = ylim, main = "SPCC")
  } else if ("SPRCC" %in% names(x)) {
    nodeplot(x$SPRCC, ylim = ylim, main = "SPRCC")
  }
}


ggplot.pcc <- function(x, ylim = c(-1,1), ...) {  
  if ("PCC" %in% names(x)) {
    nodeggplot(listx = list(x$PCC), xname = "PCC", ylim = ylim, title = "PCC")
  } else if ("PRCC" %in% names(x)) {
    nodeggplot(listx = list(x$PRCC), xname = "PRCC", ylim = ylim, title = "PRCC")
  } else if ("SPCC" %in% names(x)) {
    nodeggplot(listx = list(x$SPCC), xname = "SPCC", ylim = ylim, title = "SPCC")
  } else if ("SPRCC" %in% names(x)) {
    nodeggplot(listx = list(x$SPRCC), xname = "SPRCC", ylim = ylim, title = "SPRCC")
  }
}

