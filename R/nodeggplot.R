
#                       Nodeggplot: anti-boxplot in ggplot
#                         Sebastien Da Veiga (June 2019, Seignosse)

if (getRversion() >= "2.15.1") utils::globalVariables(c("x","y","id"))

nodeggplot <- function(listx, xname, xlim = NULL, ylim = NULL, labels = TRUE, title = NULL, 
                       col = par("col"), pch = 21, at = NULL, y_col = NULL, y_dim3 = NULL, ...) {
  
  ngraphs <- length(listx)
  x <- unlist(listx)
  n <- nrow(listx[[1]])
  if (is.null(xlim)) {
    xlim <- c(1, n)
  }
  if (n<=10){
    angle <- 0
    hjust <- 0
  }
  if (n>10 & n <=20){
    angle <- 45
    hjust <- 1
  }
  if (n>20){
    angle <- 90
    hjust <- 1
  }
  if (is.null(ylim)) {
    ylim <- c(min(x), max(x))
  }
  if (is.null(at)) {
    at <- 1 : n
  }
  if (is.null(title)){
    title <- title
  }
  
#  if (class(labels) == "logical") {
  if (inherits(labels, "logical")){
    if (labels) {
      l <- rownames(listx[[1]])
    } else {
      l <- NULL
    }
#  } else if (class(labels) == "character") {
  } else if (inherits(labels, "character")){
    l <- labels
  }
  
  # bias
  
  d <- NULL
  for (i in 1:ngraphs){
    x <- listx[[i]]
    if ("bias" %in% dimnames(x)[[2]]) {
      if(is.null(y_col) && is.null(y_dim3)){
        xx <- x[["original"]] - x[["bias"]]
      } else if(!is.null(y_col) && is.null(y_dim3)){
        xx <- x[, "original", y_col] - x[, "bias", y_col]
      } else if(!is.null(y_col) && !is.null(y_dim3)){
        xx <- x[, "original", y_col, y_dim3] - x[, "bias", y_col, y_dim3]
      }
    } else {
      if(is.null(y_col) && is.null(y_dim3)){
        xx <- x[["original"]]
      } else if(!is.null(y_col) && is.null(y_dim3)){
        xx <- x[, y_col]
      } else if(!is.null(y_col) && !is.null(y_dim3)){
        xx <- x[, y_col, y_dim3]
      }
    }
    d <- rbind(d,data.frame(x=at,y=xx,id=xname[i]))
  }
  
  # confidence intervals
  d2 <- NULL
  n2 <- 0
  for (i in 1:ngraphs){
    x <- listx[[i]]
    if (("min. c.i." %in% dimnames(x)[[2]]) & "max. c.i." %in% dimnames(x)[[2]]) {
      if(is.null(y_col) && is.null(y_dim3)){
        min_ci <- x[["min. c.i."]]
        max_ci <- x[["max. c.i."]]
      } else if(!is.null(y_col) && is.null(y_dim3)){
        min_ci <- x[, "min. c.i.", y_col]
        max_ci <- x[, "max. c.i.", y_col]
      } else if(!is.null(y_col) && !is.null(y_dim3)){
        min_ci <- x[, "min. c.i.", y_col, y_dim3]
        max_ci <- x[, "max. c.i.", y_col, y_dim3]
      }
      n2 <- n2 +1
    }else{
      min_ci <- rep(NA,n)
      max_ci <- rep(NA,n)
    }
    d2 <- rbind(d2,data.frame(min_ci=min_ci,max_ci=max_ci))
  }
  
  d <- cbind(d,d2)
  
  if (ngraphs>1){
    pd <- position_dodge(0.3)
    if (n2>0){
      g <- ggplot(d, aes(x=x, y=y, shape=id)) + 
        geom_point(size=3, colour=col, position = pd) +
        geom_errorbar(aes(ymin=min_ci, ymax=max_ci), width=.1, position = pd) + 
        scale_shape_manual(values=pch) +
        coord_cartesian(ylim=ylim) +
        labs(y="", x = "", shape = title) +
        scale_x_discrete(limits=l) + 
        theme_bw() +
        theme(axis.text.x = element_text(face = "bold", size = 12, angle = angle, hjust = hjust), axis.text.y = element_text(face = "bold", size = 12), legend.position = c(0.8, 0.9))
    }else{
      g <- ggplot(d, aes(x=x, y=y, shape=id)) + 
        geom_point(size=3, colour=col, position = pd) +
        scale_shape_manual(values=pch) +
        coord_cartesian(ylim=ylim) +
        labs(y="", x = "", shape = title) +
        scale_x_discrete(limits=l) + 
        theme_bw() +
        theme(axis.text.x = element_text(face = "bold", size = 12, angle = angle, hjust = hjust), axis.text.y = element_text(face = "bold", size = 12), legend.position = c(0.8, 0.9))
    }
  }else{
    if (n2>0){
      g <- ggplot(d, aes(x=x, y=y)) + 
        geom_point(size=3, shape=pch, colour=col) +
        geom_errorbar(aes(ymin=min_ci, ymax=max_ci), width=.1, colour=col) + 
        coord_cartesian(ylim=ylim)+
        labs(title= xname, y="", x = "") +
        scale_x_discrete(limits=l) + 
        theme_bw() +
        theme(axis.text.x = element_text(face = "bold", size = 12, angle = angle, hjust = hjust), axis.text.y = element_text(face = "bold", size = 12), plot.title = element_text(face = "bold", size = 15))
    }else{
      g <- ggplot(d, aes(x=x, y=y)) + 
        geom_point(size=3, shape=pch, colour=col) +
        coord_cartesian(ylim=ylim)+
        labs(title= xname, y="", x = "") +
        scale_x_discrete(limits=l) + 
        theme_bw() +
        theme(axis.text.x = element_text(face = "bold", size = 12, angle = angle, hjust = hjust), axis.text.y = element_text(face = "bold", size = 12), plot.title = element_text(face = "bold", size = 15)) 
    }
  }
  return(g)
}
