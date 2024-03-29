# Standardized Regression Coefficients
#
# Gilles Pujol 2006
# Bertrand Iooss 2020 for logistic model


estim.src <- function(data, logistic, i = 1:nrow(data) ) {
  d <- data[i, ]
  if (!logistic){
    lm.Y <- lm(formula(paste(colnames(d)[1], "~", paste(colnames(d)[-1], collapse = "+"))), data = d)
    coefficients(lm.Y)[-1] * sapply(d[-1], sd) / sapply(d[1], sd)
  } else{
    glm.Y <- glm(formula(paste(colnames(d)[1], "~", paste(colnames(d)[-1], collapse = "+"))), family = "binomial", data = d)
    varY <- glm.Y$linear.predictors/(1-glm.Y$deviance/glm.Y$null.deviance)
    coefficients(glm.Y)[-1] * sapply(d[-1], sd) / sd(varY)
  }
}

src <- function(X, y, rank = FALSE, logistic = FALSE, nboot = 0, conf = 0.95) {
  data <- data.frame(Y = y, X)
  
  if (logistic) rank <- FALSE # Impossible to perform logistic regression with a rank transformation
  
  if (rank) {
    for (i in 1:ncol(data)) {
      data[,i] <- rank(data[,i])
    }
  }
  
  if (nboot == 0) {
    src <- data.frame(original = estim.src(data, logistic ))
    rownames(src) <- colnames(X)
  } else {
    boot.src <- boot(data, estim.src, logistic = logistic, R = nboot)
    src <- bootstats(boot.src, conf, "basic")
    rownames(src) <- colnames(X)
  }

  out <- list(X = X, y = y, rank = rank, nboot = nboot, conf = conf,
              call = match.call())
  class(out) <- "src"
  if (! rank) {
    out$SRC <- src
  } else {
    out$SRRC = src
  }
  return(out)
}


print.src <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if ("SRRC" %in% names(x)) {
    cat("\nStandardized Rank Regression Coefficients (SRRC):\n")
    print(x$SRRC)
  } else if ("SRC" %in% names(x)) {
    cat("\nStandardized Regression Coefficients (SRC):\n")
    print(x$SRC)
  }
}


plot.src <- function(x, ylim = c(-1,1), ...) {  
  if ("SRRC" %in% names(x)) {
    nodeplot(x$SRRC, ylim = ylim, main = "SRRC")
  } else if ("SRC" %in% names(x)) {
    nodeplot(x$SRC, ylim = ylim, main = "SRC")
  }
}

ggplot.src <- function(data, mapping = aes(), ylim = c(-1,1), ..., environment = parent.frame()) {  
  x <- data
  if ("SRRC" %in% names(x)) {
    nodeggplot(listx = list(x$SRRC), xname = "SRRC", ylim = ylim, title = "SRRC")
  } else if ("SRC" %in% names(x)) {
    nodeggplot(listx = list(x$SRC), xname = "SRC", ylim = ylim, title = "SRC")
  }
}
