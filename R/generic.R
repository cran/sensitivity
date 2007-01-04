
## File: generic.R
## Description: functions that give some genericity
## Author: Gilles Pujol


tell <- function(sa, y = NULL)
  UseMethod("tell")


response <- function(sa, ...){
  id <- deparse(substitute(sa))
  
  if (class(sa$model) == "function")
     y <- sa$model(sa$x, ...)
   else if (TRUE %in% (paste("predict.", class(sa$model), sep="") %in%
                       methods(predict)))
     y <- predict(sa$model, sa$x, ...)
   else
     stop("The model isn't a function or a predictor")

  if (class(y) != "numeric"){
    y <- as.numeric(y)
    warning("Conversion of the response to numeric")
  }

  sa$y <- y
  assign(id, sa, parent.frame())
}


sobol <- function(method = "sobol93", model = NULL, x1, x2, max.order = 1,
                  nboot = 0, conf = 0.95, ...){
  switch(method,
         sobol93 = sobol.sobol93(model, x1, x2, max.order, nboot, conf, ...),
         saltelli02 = sobol.saltelli02(model, x1, x2, nboot, conf, ...)
         ) 
}


fast <- function(method = "saltelli99", model = NULL, factors, n, M = 4,
                 omega = NULL, q = NULL, q.arg = NULL, ...){
  switch(method,
         saltelli99 = fast.saltelli99(model, factors, n, M, omega, q, q.arg, ...)
         )
}
