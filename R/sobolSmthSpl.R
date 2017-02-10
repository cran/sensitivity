# Author: Filippo Monari

sobolSmthSpl <- function (Y, X) {
  #Determines the Si coefficient for singular parameters through smoothing with roughness penalty.
  #Reference: Saltelli, A; Ratto, M; Andres, T; Campolongo, F; Cariboni, J; Gatelli, D; Saisana, M & Tarantola, S.
  #Global Sensitivity Analysis: The Primer Wiley-Interscience, 2008
  #Arguments:
  #Y: matrix of model outputs (only one column)
  #X: matrix model parameters
  #MAIN
  ANS = list()
  ANS[['call']] = match.call()
  ANS[['X']] = X
  ANS[['Y']] = Y
  par.names = colnames(X)	#gets parameters names
  if (is.null(colnames(X))) par.names = paste0('X', 1:ncol(X))
  X = normalize(X)    		#normalize inpiuts between [0,1]
  Y = Y - mean(Y)     		#center model responses
  Y = sapply(1:ncol(X), function(i)return(Y[order(X[,i])]))    #order Y before X (or create a new variable for X)
  X = sapply(1:ncol(X), function(i)return(X[order(X[,i]),i]))   
  SMTH = optSmooth(Y, X, c(-2, 2))
  SA.tab = t(sapply(SMTH, est.Si))
  colnames(SA.tab) = c('Si', 'se', 'q0.05')
  rownames(SA.tab) = par.names
  ANS[['S']] = SA.tab 
  class(ANS) = 'sobolSmthSpl'
  return(ANS)
}

est.Si <-function (SMTH) {
  gi = SMTH[['y']]
  yi = SMTH[['yin']]
  Si = var(gi) / var(yi)
  #calculates the standard error of the main effect estimates
  yi.sc = (yi - mean(yi)) / sd(yi)            	#scaled yi
  u = (yi - gi) / (sd(yi) * abs(1 - Si)**0.5)   #scales residuals
  Si.se = abs(1 - Si) * sd(yi.sc**2 - u**2) / length(yi)**0.5
  q0.05 = qnorm(0.05, Si, Si.se)
  return(c(Si, Si.se, q0.05))
}

optSmooth <- function (Y, X, interval) {
  #DOC
  #Optimises the spline smoothing across all the column of the matrix Y.
  #Argumnts:
  #Y: matrix containing the observation to smooth
  #X: matrix of inputs
  #interval: interval where to search for the spar smoothing parameter
  #CONTAINS
  objective = function(spar) {
    nc = min(ncol(Y), parallel::detectCores())
    ANS = try(parallel::mclapply(1:ncol(Y), function(i)smooth.spline(X[,i], Y[,i], spar = spar, all.knots = T), mc.cores = nc), silent = T)
    if (is(ANS, 'try-error')) { #Windows does not support mclapply...
		SMTH <<- lapply(1:ncol(Y), function(i)smooth.spline(X[,i], Y[,i], spar = spar, all.knots = T))
	} else {
		SMTH <<- ANS
	}
    return(sum(sapply(SMTH, function(x)x$cv)))
  }
  SMTH = NULL
  ANS = optimize(f = objective, interval = interval, tol = sqrt(.Machine$double.eps), maximum = FALSE)
  return(SMTH)
}

normalize <- function (X, MAXS = NULL, MINS = NULL, inv = F) {
	#DOC
	#Normalizes a vector or a matrix according to the given 'MAXS' and 'MINS'.
	#If 'MAXS' and 'MINS' are not provided 'X' is scaled so that each column is within [0,1].
	#ARGUMENTS
	#X: matrix or vector to normalize
	#MAXS, MINS: maxima and minima. NULL or NA values are set to max(X[,i]) and min(X[,i]) respectively.
	#inv: if T performs the inverse transformation
	#MAIN
	X = as.matrix(X)
	if (inv) {
		return(scale(scale(X, center = F, scale = (MAXS - MINS)**(-1)), center = -MINS, scale = F))
	} else {
		if (is.null(MAXS)) MAXS = apply(X, 2, max)
		if (is.null(MINS)) MINS = apply(X, 2, min)
		for (i in 1:ncol(X)) {
			if (is.na(MAXS[i])) MAXS[i] = max(X[,i])
			if (is.na(MINS[i])) MINS[i] = min(X[,i])
		}
		return(scale(X, center = MINS, scale = MAXS - MINS))
	}
}	

plot.sobolSmthSpl <- function(x, ...) {
	yrng = range(c(1 + x[['S']][,'se'], x[['S']][,'q0.05']))
	#plot estimates
	plot(x = 1:nrow(x[['S']]), y = x[['S']][,'Si'], pch = 19, ylim = yrng, xaxt = 'n', xlab = 'parameter', ylab = 'Si and se', ...)
	#plot q 0.05
	points(x = 1:nrow(x[['S']]), y = x[['S']][,'q0.05'], pch = 19, col = 2)
	#plot se
	arrows(x0 = 1:nrow(x[['S']]), y0 = x[['S']][,'Si'], y1 = x[['S']][,'Si'] - x[['S']][,'se'], length = 0)	#lower
	arrows(x0 = 1:nrow(x[['S']]), y0 = x[['S']][,'Si'], y1 = x[['S']][,'Si'] + x[['S']][,'se'], length = 0)	#upper
	#plot 0
	abline(h = 0, lty = 2)
	#x axis
	axis(side = 1, at = 1:nrow(x[['S']]), labels = row.names(x[['S']]))
	#legend
	legend('topright', legend = c('Si', 'se', 'q0.05'), pch = c(19, NA, 19), lty = c(NA, 1, NA), col = c(1, 1, 2), horiz = T, bty = 'n')
}

print.sobolSmthSpl <- function(x, ...) {
	cat("\nCall:\n", deparse(x[['call']]), "\n", sep = "")
	cat("\nModel runs:", length(x[['Y']]), "\n")
    cat("\nFirst order indices:\n")
    print(x[['S']])
}	
