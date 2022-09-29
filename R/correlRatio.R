# Correlation ratio between a quantitative variable and a qualitative one
#
# Bertrand Iooss 2022 (inspired from fct corRatio() of the DiscriMiner package)


correlRatio <-
  function(X, y)
  {
    # Correlation ratio
    # X: sample of a quantitative variable
    # y: sample of factor variables
    
    if (!is.numeric(X)) 
      stop("\n'X' must be a numeric vector")
    if (!is.factor(y)) y = as.factor(y)
    if (nlevels(y) == 1)
      stop("\n'y' has only one category")
    # correlation ratio
    regr <- lm(X ~ y)
    res <- summary(regr)$r.squared
    res
  }