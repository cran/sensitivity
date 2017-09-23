# Author : Olivier Roustant (2017)

support <- function(model, X, Xnew = NULL, fX = NULL, gradfX = NULL, h = 1e-6, ...) {

  # column j of gradfX contains (a finite difference approx. of)
  # the derivative with respect to Xj evaluated at the sample matrix X

  n <- nrow(X)
  if (is.null(Xnew)) Xnew <- X
  n.points <- nrow(Xnew)
  d <- ncol(X)

  main <- total <- matrix(NA, n.points, d)

    # evaluate f on X
  if (is.null(fX)) fX <- model(X, ...)
  if (is.null(gradfX)) {
    gradfX <- matrix(NA, n, d)
    missingGrad <- TRUE
  }

  for (j in 1:d){
      if (missingGrad){
        # compute finite differences
        Xplus <- X
        Xplus[, j] <- Xplus[, j] + h
        gradfX[,j] <- (model(Xplus) - fX)/h
      }

    DGSM <- colMeans(gradfX^2)
    # compute non-parametric estimate of conditional expectations
    mydata <- data.frame(x = X[, j], y = gradfX[, j], y2 = gradfX[, j]^2)
    newdata <- Xnew[, j]
    m1 <- smooth.spline(x = mydata$x, y = mydata$y)
    m2 <- smooth.spline(x = mydata$x, y = mydata$y2)
    main[, j] <- (predict(m1, newdata)$y)^2
    total[, j] <- predict(m2, newdata)$y
  }

  x <- list(main = main, total = total, DGSM = DGSM, 
            X = X, Xnew = Xnew, fX = fX, gradfX = gradfX)
  class(x) <- "support"
  return(x)
}


arg2arg <- function(p, p.arg, d){
  if (length(p) == 1) p <- rep(p, d)
  if (is.null(p.arg)) {
    p.arg <- rep(list(list()), d)
  } else if (FALSE %in% sapply(p.arg, is.list)) {
    p.arg <- rep(list(p.arg), d)
  }
  list(p = p, p.arg = p.arg)
}


plot.support <- function(x, i = 1:ncol(x$X),
                         xprob = FALSE, p = NULL, p.arg = NULL,
                         ylim = NULL, col = 1:3, lty = 1:3, lwd=c(2,2,1), cex = 1, ...){

  d <- ncol(x$main)
  S <- sum(x$DGSM)
  indices <- i
  
  if ((xprob) & (!is.null(p))){
    res <- arg2arg(p, p.arg, d)
    p <- res$p
    p.arg <- res$p.arg
  }

  if (is.null(ylim)) ylim = c(0, max(x$total)/S)

  for (j in indices){
    xj <- x$X[, j]
    xjnew <- x$Xnew[, j]
    if (xprob) {   # transform the x-axis to probabilities P(Xj < x)
      if (is.null(p)){  # if p is not given use the empirical cdf
        xjnew <- ecdf(xj)(xjnew)
#        xj <- ecdf(xj)(xj)
      } else {
        xjnew <- do.call(p[j], c(list(xjnew), p.arg[[j]]))
#        xj <- do.call(p[j], c(list(xj), p.arg[[j]]))
      }
    }
    ind <- sort(xjnew, index.return = TRUE)$ix
    plot(xjnew[ind], x$total[ind, j]/S, type = "l", lty = lty[1],
         xlab = substitute(x[j], list(j = j)),
         ylab="",
         col = col[1], ylim = ylim, lwd = lwd[1])
    if (xprob) title(sub = '(probability scale)')
#    points(xj, x$gradfX[,j]^2/S, pch = 19, cex = 0.2, col = "grey")
    lines(xjnew[ind], x$main[ind, j]/S, col = col[2], lty = lty[2], lwd = lwd[2])
    abline(h = x$DGSM[j]/S, lty = lty[3], col = col[3], lwd = lwd[3])
    legend('topleft', c(expression(zeta ** symbol("*T*")),
                        expression(zeta ** symbol("*")),
                        expression(nu ** symbol("*"))),
           col = col, lty = lty, lwd = lwd, bty = "n",
           horiz = FALSE, adj = 0, cex = cex)
  
}

}


scatterplot.support <- function(x, i = 1:ncol(x$X), 
                                xprob = FALSE, p = NULL, p.arg = NULL, 
                                cex = 1, cex.lab = 1, ...){

  d <- ncol(x$main)
  S <- sum(x$DGSM)
  indices <- i
  
  if ((xprob) & (!is.null(p))){
    res <- arg2arg(p, p.arg, d)
    p <- res$p
    p.arg <- res$p.arg
  }
  
  for (j in indices){
    xj <- x$X[, j]
    if (xprob) {   # transform the x-axis to probabilities P(Xj < x)
      if (is.null(p)){  # if p is not given use the empirical cdf
        xj <- ecdf(xj)(xj)
      } else {
        xj <- do.call(p[j], c(list(xj), p.arg[[j]]))
      }
    }
    par(mar = c(5, 5, 4, 2) + 0.1)   # enlarge the left margin of 1 unit
    df <- data.frame(x = xj, y = x$gradfX[, j]/S)
    if (requireNamespace("ggplot2", quietly = TRUE)){
      p1 <- ggplot2::ggplot(df, ggplot2::aes(x, y)) + ggplot2::geom_point(size = 0.5) +
      ggplot2::theme_bw() + #removeGrid()+
      ggplot2::labs(x = substitute(x[j], list(j = j)), y = substitute(df/dx[j], list(j = j))) +
      ggplot2::theme(axis.title = ggplot2::element_text(size = 10*cex.lab)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10*cex)) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10*cex))
    if (xprob) {
      p1 <- p1 + ggplot2::labs(caption = '(probability scale)')
      print(p1)
    } else {
      if (requireNamespace("ggExtra", quietly = TRUE)){ 
      p2 <- ggExtra::ggMarginal(p1, margins = "x", ...)
      print(p2)}
    }}
  }
}

