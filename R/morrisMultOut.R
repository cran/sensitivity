
morrisMultOut <- function (model = NULL, factors, r = 50, design = list(type = "oat", levels = 5, grid.jump = 3), binf = 0, bsup = 1, scale = TRUE, ...) {
	M = morris(model = NULL, factors = factors, r = r, design = design, binf = binf, bsup = bsup, scale = scale, ...)
	if (!is.null(model)) {
		Y = model(M$X)
		M = .morrisMultOut(Y, M)
	}
	class(M) = c('morrisMultOut', 'morris')
	return(M)
}

.morrisMultOut <- function (Y, M, ...) {
#	class(M) = 'morris'	#delete
	SVD = svd(Y)
	W = SVD$d**2 / sum(SVD$d**2)
	ee = 0
	for (i in 1:ncol(SVD$v)) {
		tell.morris(M, SVD$v[,i])	
		ee = ee + M$ee**2 * W[i]
	}
	M$ee = sqrt(ee)
#	class(M) = c('morrisMultOut', 'morris')	
	return(M)
}

tell.morrisMultOut <- function(x, y = NULL, ...) {
	id <- deparse(substitute(x))	
	ANS = .morrisMultOut(y, x)
	assign(id, ANS, parent.frame())
}


