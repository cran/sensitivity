# Kriging-based global sensitivity analysis taking into account both 
# the meta-model and the Monte-Carlo errors

# Author : loic le Gratiet, 2014

sobolGP <- function (
model,  
type="SK",
MCmethod="sobol",                                                                                                                                                               
X1,  
X2, 
nsim=100, 
nboot=1,
conf = 0.95,
sequential = FALSE, 
candidate = NULL, 
sequential.tot=FALSE,
max_iter = 1000
) 
{
    if(sequential){
    	ncandidate <- dim(candidate)[1]
    	dcandidate <- dim(candidate)[2]
    }	

    if ((ncol(X1) != ncol(X2)) | (nrow(X1) != nrow(X2))) 
        stop("The samples X1 and X2 must have the same dimensions")
    if(sequential){
	if(is.null(candidate)){
	  		stop("The number of candidate points must be greater than zero")
	}
    	if (ncol(X1) != ncol(candidate)){
	  		stop("The candidate points, X1 and X2 must have the same dimensions")
	}
    }
    if(MCmethod=="sobol"||MCmethod=="sobolEff"){
	if(sequential.tot){
	  stop("Sequential design for total indices is only available for sobol2002, sobol2007 and soboljansen methods")
	}
    }

    p <- ncol(X1)

    S <- list()

    if(sequential){
	Svar <- matrix(nrow = ncandidate, ncol = p)
    }
    if(sequential.tot){
	STotvar <- matrix(nrow = ncandidate, ncol = p)
    }

    output <- list()
    output$call$X1 <- X1
    output$call$X2 <- X2
    output$call$conf <- conf
    output$call$nboot <- nboot
    output$call$candidate <- candidate
    output$call$sequential <- sequential
    output$call$max_iter <- max_iter
    output$call$sequential.tot <- sequential.tot
    output$call$model <- model

Tot=FALSE
if(MCmethod=="sobol2002"||MCmethod=="sobol2007"||MCmethod=="soboljansen") Tot=TRUE

if(MCmethod!="sobol2007"){

    for (i in 1:p) {
        Xb <- X2
        Xb[, i] <- X1[, i]
        X <- rbind(X1, Xb)
	  nX <- dim(X)[1]

	if(sequential){
	  X <- rbind(X, data.frame(candidate))
	}
		
	rm(list=c("Xb"))

	  ysimu <- simulateGP.sobol(object = model, nsim = nsim,  newdata=X, 
                            cond=TRUE, checkNames=FALSE, max_iter=1000,type)

	if(MCmethod=="sobol"||MCmethod=="sobol2002"){
	  	S[[i]] <- sobolpickfreeze(ysimu[,1:(nX/2)] , ysimu[,(nX/2+1):nX],nboot)
	}
	if(MCmethod=="sobolEff"){
		S[[i]] <- sobolEffpickfreeze(ysimu[,1:(nX/2)] , ysimu[,(nX/2+1):nX],nboot)
	}
	if(MCmethod=="soboljansen"){
		S[[i]] <- soboljansenpickfreeze(ysimu[,1:(nX/2)] , ysimu[,(nX/2+1):nX],nboot)
	}

	  if(sequential){
	  	predCov <- predictGP.sobol(object = model, newdata1=data.frame(X), newdata2=data.frame(candidate), type=type, prednewdata1 = FALSE, prednewdata2 = TRUE)

	  	for(k in 1:ncandidate){
			ynew <- predCov$mean2[k]
			zsimu <- t(as.matrix(predCov$cov[-c((nX+1):(nX+ncandidate)),k])%*%rep(1,nsim))*(ynew-ysimu[,(nX+k)])/predCov$cov[(nX+k),k]+ysimu[,-c((nX+1):(nX+ncandidate))]

			if(MCmethod=="sobol"||MCmethod=="sobol2002"){
	  			Scand <- sobolpickfreeze(zsimu[,1:(nX/2)] , zsimu[,(nX/2+1):nX],nboot=1) 
			}
			if(MCmethod=="sobolEff"){
				Scand <- sobolEffpickfreeze(zsimu[,1:(nX/2)] , zsimu[,(nX/2+1):nX],nboot=1) 
			}
			if(MCmethod=="soboljansen"){
				Scand <- soboljansenpickfreeze(zsimu[,1:(nX/2)] , zsimu[,(nX/2+1):nX],nboot=1) 
			}

			rm(list=c("zsimu"))
			Svar[k,i] <- var(Scand)
			rm(list=c("Scand"))
	  	}
		rm(list=c("predCov", "ynew"))
	  }
	rm(list=c("ysimu" ))
    }

	

	if(Tot){
	
	    STot <- list()
	
    	for (i in 1:p) {
	        Xb <- X1
	        Xb[, i] <- X2[, i]
        	X <- rbind(X1, Xb)
		  nX <- dim(X)[1]
	  	if(sequential.tot){
		  X <- rbind(X, data.frame(candidate))
	  	}
			
		rm(list=c("Xb"))
	
		  ysimu <- simulateGP.sobol(object = model, nsim = nsim,  newdata=X, 
	                            cond=TRUE, checkNames=FALSE, max_iter=1000,type)
	
		if(MCmethod=="sobol2002"){ 
			STot[[i]] <- sobolT2002pickfreeze(ysimu[,1:(nX/2)],ysimu[,(nX/2+1):nX],nboot)
		}
		if(MCmethod=="soboljansen"){ 
			STot[[i]] <- sobolTjansenpickfreeze(ysimu[,1:(nX/2)],ysimu[,(nX/2+1):nX],nboot)
		}
	
		  if(sequential.tot){

		  	predCov <- predictGP.sobol(object = model, newdata1=data.frame(X), newdata2=data.frame(candidate), type=type)
	
		  	for(k in 1:ncandidate){
				ynew <- predCov$mean[k]
				zsimu <- t(as.matrix(predCov$cov[-c((nX+1):(nX+ncandidate)),k])%*%rep(1,nsim))*(ynew-ysimu[,(nX+k)])/predCov$cov[(nX+k),k]+ysimu[,-c((nX+1):(nX+ncandidate))]

				if(MCmethod=="sobol2002"){ 
					STotcand <- sobolT2002pickfreeze(zsimu[,1:(nX/2)],zsimu[,(nX/2+1):nX],nboot=1)
				}
				if(MCmethod=="soboljansen"){ 
					STotcand <- sobolTjansenpickfreeze(zsimu[,1:(nX/2)],zsimu[,(nX/2+1):nX],nboot=1)
				}

				rm(list=c("zsimu"))
				STotvar[k,i] <- var(STotcand)
				rm(list=c("STotcand"))
		  	}
		  rm(list=c("predCov", "ynew"))
	  	}
		rm(list=c("ysimu" ))
	    }
	}

		rm(list=c("X", "X1","X2"))
} 

if(MCmethod=="sobol2007"){
    STot <- list()

    for (i in 1:p) {
        Xb <- X1
        Xb[, i] <- X2[, i]
        X <- rbind(X1, Xb, X2)
	nX <- dim(X)[1]

	if(sequential||sequential.tot){
		X <- rbind(X, data.frame(candidate))
	}		

	rm(list=c("Xb"))

	  ysimu <- simulateGP.sobol(object = model, nsim = nsim,  newdata=X, 
                            cond=TRUE, checkNames=FALSE, max_iter=1000,type)
	

	S[[i]] <- sobol2007pickfreeze(ysimu[,1:(nX/3)] , ysimu[,(nX/3+1):(2*nX/3)] , ysimu[,(2*nX/3+1):nX],nboot)
	STot[[i]] <- sobolT2007pickfreeze(ysimu[,1:(nX/3)] ,ysimu[,(nX/3+1):(2*nX/3)] ,nboot)	

	  if(sequential||sequential.tot){

	  	predCov <- predictGP.sobol(object = model, newdata1=data.frame(X), newdata2=data.frame(candidate), type=type)

	  	for(k in 1:ncandidate){
			ynew <- predCov$mean[k]
			zsimu <- t(as.matrix(predCov$cov[-c((nX+1):(nX+ncandidate)),k])%*%rep(1,nsim))*(ynew-ysimu[,(nX+k)])/predCov$cov[(nX+k),k]+ysimu[,-c((nX+1):(nX+ncandidate))]

			if(sequential){
	  			Scand <- sobol2007pickfreeze(zsimu[,1:(nX/3)] , zsimu[,(nX/3+1):(2*nX/3)] , zsimu[,(2*nX/3+1):nX],nboot=1)
				Svar[k,i] <- var(Scand)
				rm(list=c("Scand"))
			}
			if(sequential.tot){
				STotcand <- sobolT2007pickfreeze(zsimu[,1:(nX/3)] ,zsimu[,(nX/3+1):(2*nX/3)], nboot=1 )	
				STotvar[k,i] <- var(STotcand)
				rm(list=c("STotcand"))
			}
			rm(list=c("zsimu"))

	  	}
		rm(list=c("predCov", "ynew"))
	  }
	rm(list=c("ysimu"))
    }
}

namesS <- c()
for (i in 1:p){
	namesS <- c(namesS,paste("S",i,sep=""))
}
names(S) <- namesS

if(Tot){
	namesStot <- c()
	for (i in 1:p){
		namesStot <- c(namesStot,paste("T",i,sep=""))
	}
	names(STot) <- namesStot
}


    output$S <- S
	rm(list=c("S"))

    if(sequential){
	   SumVar <- apply(Svar,1,sum)
	   output$S$xnew <- candidate[which.min(SumVar),]
	   output$S$xnewi <- which.min(SumVar)
		rm(list=c("Svar","SumVar"))
    }
   
    if(Tot){
	   output$T <- STot
		rm(list=c("STot"))
	   if(sequential.tot){
	   	SumVar <- apply(STotvar,1,sum)
	   	output$T$xnew <- candidate[which.min(SumVar),]
	   	output$T$xnewi <- which.min(SumVar)	
		rm(list=c("STotvar","SumVar"))	
	   }
    }

	class(output) <- "sobolGP"

	output$S$mean <- matrix(nrow = 1, ncol = p)
	output$S$var <- matrix(nrow = 1, ncol = p)
	output$S$ci <- matrix(nrow = 2, ncol = p)

	output$S$varPG <- matrix(nrow = 1, ncol = p)
	output$S$varMC <- matrix(nrow = 1, ncol = p)		

	if(Tot){
		output$T$mean <- matrix(nrow = 1, ncol = p)
		output$T$var <- matrix(nrow = 1, ncol = p)
		output$T$ci <- matrix(nrow = 2, ncol = p)

		output$T$varPG <- matrix(nrow = 1, ncol = p)
		output$T$varMC <- matrix(nrow = 1, ncol = p)	
	}
	for(i in 1:p){
		output$S$mean[1,i] <- mean(as.numeric(output$S[[i]]))
		output$S$var[1,i] <- var(as.numeric(output$S[[i]]))
		output$S$ci[1,i] <- quantile(as.numeric(output$S[[i]]), (1-conf)/2) 
		output$S$ci[2,i] <- quantile(as.numeric(output$S[[i]]), (1+conf)/2) 

		if(nboot==1){
			output$S$varPG[1,i] <- var(output$S[[i]])
		} else {
			output$S$varPG[1,i] <- mean(apply(output$S[[i]],1,var))
			output$S$varMC[1,i] <- mean(apply(output$S[[i]],2,var))
		}

		if(Tot){
			output$T$mean[1,i] <- mean(as.numeric(output$T[[i]]))
			output$T$var[1,i] <- var(as.numeric(output$T[[i]]))
			output$T$ci[1,i] <- quantile(as.numeric(output$T[[i]]), (1-conf)/2) 
			output$T$ci[2,i] <- quantile(as.numeric(output$T[[i]]), (1+conf)/2) 

			if(nboot==1){
				output$T$varPG[1,i] <- var(output$T[[i]])
			} else {
				output$T$varPG[1,i] <- mean(apply(output$T[[i]],1,var))
				output$T$varMC[1,i] <- mean(apply(output$T[[i]],2,var))
			}
		}
	}

	output$call$tot <- Tot
	output$call$method <- MCmethod
	output$call$type <- type
	output$call$nsim <- nsim
	

	return(output)
}

###############################################################################

plot.sobolGP <- function(x, ...){
  sGP <- x
  dev.new(width = 12, height = 5)
  if(sGP$call$tot){
    par(mfrow=c(1,3))
  } else {
    par(mfrow=c(1,2))
  }
  d <- length(sGP$S$mean)
  ylim=c(-0.05,1.05)
  xlim=c(1,d+0.1)
  plot(1:d,sGP$S$mean, xlim = xlim, ylim = ylim, axes = FALSE,
       xlab = "", ylab = "", pch = 21)
  axis(side = 1, at = 1:d, labels = 1:d)
  axis(side = 2, at = seq(0,1,0.1), labels = seq(0,1,0.1))
  segments(1:d,sGP$S$ci[1,],1:d,sGP$S$ci[2,])
  box()
  if(sGP$call$tot){
    points((1:d)+0.1,sGP$T$mean, pch = 24)	
    legend("topright", legend = c("main effect", "total effect"), pch = c(21,24))
    segments(1:d+0.1,sGP$T$ci[1,],1:d+0.1,sGP$T$ci[2,],lty=2)
  } else {
    legend("topright", legend = c("main effect"), pch = 21)
  }
  abline(0,0)
  
  nboot <- sGP$call$nboot
  if(nboot!=1){
    data <- data.frame(rbind(sGP$S$varPG,sGP$S$varMC))
    names(data) <- names(sGP$S)[1:d]
    barplot(as.matrix(data), angle = c(45,-45), density = 20, col = "black",
            legend = c("GP variance","MC variance"),ylim=c(0,1.2*max(sGP$S$varPG+sGP$S$varMC)))
    title(main = list("Variance decomposition of the main effects", font = 4))
    
    if(sGP$call$tot){
      data <- data.frame(rbind(sGP$T$varPG,sGP$T$varMC))
      names(data) <- names(sGP$T)[1:d]
      barplot(as.matrix(data), angle = c(45,-45), density = 20, col = "black",
              legend = c("GP variance","MC variance"),ylim=c(0,1.2*max(sGP$T$varPG+sGP$T$varMC)))
      title(main = list("Variance decomposition of the total effects", font = 4))
    }
  } else {
    data <- data.frame(sGP$S$varPG)
    names(data) <- names(sGP$S)[1:d]
    barplot(as.matrix(data), angle = 45, density = 20, col = "black",
            legend = c("GP variance"),ylim=c(0,1.2*max(sGP$S$varPG)))
    title(main = list("Variance of the main effects", font = 4))
    
    if(sGP$call$tot){
      data <- data.frame(sGP$T$varPG)
      names(data) <- names(sGP$T)[1:d]
      barplot(as.matrix(data), angle = 45, density = 20, col = "black",
              legend = c("GP variance"),ylim=c(0,1.2*max(sGP$T$varPG)))
      title(main = list("Variance of the total effects", font = 4))
    }
  }
  
}  

######################################

print.sobolGP <- function(x, ...){
  sGP <- x
  cat("\nMethod: ", sGP$call$method, "\n", sep = "")
  cat("\nModel runs:", sGP$call$model@n, "\n")
  cat("\nNumber of GP realizations:", sGP$call$nsim, "\n")
  cat("\nKriging type:", sGP$call$type, "\n")
  
  d <- length(sGP$S$mean)
  df=data.frame(pointEstimate=t(sGP$S$mean), stdErr=t(sqrt(sGP$S$var)), minCI=as.matrix(sGP$S$ci[1,]), maxCI=as.matrix(sGP$S$ci[2,]))
  colnames(df)=c("estimate","std. error","min. c.i.","max. c.i.")
  rownames(df)=names(sGP$S)[1:d]
  cat("\n")
  print(df)
  cat("\n")
  if(sGP$call$tot){
    df=data.frame(pointEstimate=t(sGP$T$mean), stdErr=t(sqrt(sGP$T$var)), minCI=as.matrix(sGP$T$ci[1,]), maxCI=as.matrix(sGP$T$ci[2,]))
    colnames(df)=c("estimate","std. error","min. c.i.","max. c.i.")
    rownames(df)=names(sGP$T)[1:d]
    print(df)
    cat("\n")
  }
}

######################################

ask.sobolGP <- function(x, tot=FALSE, ...){
  sGP <- x
  output <- c()
  if(is.null(sGP$S$xnew)){
    stop("The function sobolGP must be run with the argument sequential = TRUE")
  }
  if(!tot){
    output <- sGP$S$xnew
  } else {
    if(is.null(sGP$T$xnew)){
      stop("The function sobolGP must be run with the argument sequential.tot = TRUE")
    }	
    output <- sGP$T$xnew
  }
  d <- dim(sGP$call$candidate)[2]
  return(output)
}

######################################

tell.sobolGP <- function(x, y=NULL, xpoint=NULL, newcandidate=NULL, ...){
  sGP <- x
  xx <- xpoint
  if(is.null(newcandidate)){
    d <- dim(sGP$call$candidate)[2]
    test = c(sGP$call$candidate[,1] == xx[,1])
    for(i in 2:d){
      test <- test*(sGP$call$candidate[,i] == xx[,i])
    }
    newcandidate <- sGP$call$candidate[-which(test == 1),]
  }
  
  xx <- as.matrix(xx)
  if(dim(xx)[2]==1){
    xx <- t(xx)
  }
  rownames(xx) = NULL
  m  <- update(object = sGP$call$model, newX = xx, newy = y, cov.reestim = FALSE, trend.reestim = FALSE)
  
  
  
  res <- sobolGP(
    model = m,
    type = sGP$call$type,
    MCmethod = sGP$call$method,
    X1 = sGP$call$X1,
    X2 = sGP$call$X2,
    nsim = sGP$call$nsim,
    conf = sGP$call$conf,
    nboot = sGP$call$nboot,
    sequential = sGP$call$sequential,
    candidate = newcandidate,
    sequential.tot = sGP$call$sequential.tot,
    max_iter = sGP$call$max_iter) 
  
  return(res)
  
}

######################################

predictGP.sobol <- function (
  object,
  newdata1 ,
  newdata2,
  type,
  bias.correct = FALSE,
  checkNames = FALSE, 
  prednewdata1 = FALSE, 
  prednewdata2 = TRUE) 
{
  nugget.flag <- object@covariance@nugget.flag
  X <- object@X
  y <- object@y
  if (checkNames) {
    newdata1 <- checkNames(X1 = X, X2 = newdata1, X1.name = "the design", 
                           X2.name = "newdata")
    newdata2 <- checkNames(X1 = X, X2 = newdata2, X1.name = "the design", 
                           X2.name = "newdata")
  } else {
    newdata1 <- as.matrix(newdata1)
    newdata2 <- as.matrix(newdata2)
    d.newdata <- ncol(newdata1)
    if (!identical(d.newdata, object@d)) {
      stop("newdata must have the same numbers of columns than the experimental design")
    }
    if (!identical(colnames(newdata1), colnames(X))) {
      colnames(newdata1) <- colnames(X)
    }
    if (!identical(colnames(newdata2), colnames(X))) {
      colnames(newdata2) <- colnames(X)
    }
  }
  T <- object@T
  z <- object@z
  M <- object@M
  beta <- object@trend.coef
  
  F.newdata1 <- model.matrix(object@trend.formula, data = data.frame(newdata1))
  F.newdata2 <- model.matrix(object@trend.formula, data = data.frame(newdata2))
  
  if(prednewdata1){
    y.predict.trend1 <- F.newdata1 %*% beta
  } 
  if(prednewdata2){
    y.predict.trend2 <- F.newdata2 %*% beta
  }
  if (requireNamespace("DiceKriging", quietly = TRUE)){ 
    c.newdata1 <- DiceKriging::covMat1Mat2(object@covariance, X1 = X, X2 = newdata1, 
                                           nugget.flag = object@covariance@nugget.flag)
    c.newdata2 <- DiceKriging::covMat1Mat2(object@covariance, X1 = X, X2 = newdata2, 
                                           nugget.flag = object@covariance@nugget.flag)
  }
  Tinv.c.newdata1 <- backsolve(t(T), c.newdata1, upper.tri = FALSE)
  Tinv.c.newdata2 <- backsolve(t(T), c.newdata2, upper.tri = FALSE)
  
  output.list <- list()
  if(prednewdata1){
    y.predict.complement1 <- t(Tinv.c.newdata1) %*% z
    y.predict1 <- y.predict.trend1 + y.predict.complement1
    rm(list=c("y.predict.complement1","y.predict.trend1"))
    y.predict1 <- as.numeric(y.predict1)
    output.list$mean1 <- y.predict1	
    rm(list=c("y.predict1"))
  }
  if(prednewdata2){
    y.predict.complement2 <- t(Tinv.c.newdata2) %*% z
    y.predict2 <- y.predict.trend2 + y.predict.complement2
    rm(list=c("y.predict.complement2","y.predict.trend2"))
    y.predict2 <- as.numeric(y.predict2)
    output.list$mean2 <- y.predict2
    rm(list=c("y.predict2"))
  }
  
  if (requireNamespace("DiceKriging", quietly = TRUE)){
    C.newdata <- DiceKriging::covMat1Mat2(object@covariance, X1 = newdata1, X2 = newdata2, 
                                          nugget.flag = object@covariance@nugget.flag)
  }
  
  rm(list=c("newdata1","newdata2"))
  
  cond.cov <- C.newdata - crossprod(Tinv.c.newdata1,Tinv.c.newdata2 )
  
  if (type == "UK") {
    T.M <- chol(t(M) %*% M)
    
    s2.predict.mat1 <- backsolve(t(T.M), t(F.newdata1 - 
                                             t(Tinv.c.newdata1) %*% M), upper.tri = FALSE)
    s2.predict.mat2 <- backsolve(t(T.M), t(F.newdata2 - 
                                             t(Tinv.c.newdata2) %*% M), upper.tri = FALSE)
    
    cond.cov <- cond.cov + crossprod(s2.predict.mat1,s2.predict.mat2)
    if (bias.correct) 
      cond.cov <- cond.cov * object@n/(object@n - object@p)
  }
  output.list$cov <- cond.cov
  
  return(output.list)
}

######################################

simulateGP.sobol <- function(
  object,  
  nsim = 100, 
  newdata = NULL, 
  cond = TRUE, 
  checkNames = FALSE, 
  max_iter = 1000,
  type
) {
  
  if (!is.logical(cond)) stop("'cond' must be TRUE/FALSE")
  if ((!is.null(newdata)) && (checkNames)) {
    newdata <- checkNames(X1 = object@X, X2 = newdata, X1.name = "the design", X2.name = "newdata")
  }
  
  if (is.null(newdata)) {
    newdata <- object@X
    F.newdata <- object@F
  } else {
    newdata <- as.matrix(newdata)
    m <- nrow(newdata)
    
    if (!identical(ncol(newdata), object@d)) 
      stop("newdata must have the same numbers of columns than the experimental design")
    if (!identical(colnames(newdata), colnames(object@X))) {
      colnames(newdata) <- colnames(object@X)
    }
    
    newdata <- rbind(newdata, object@X)
    row.names(newdata) <- NULL
    F.newdata <- model.matrix(object@trend.formula, data = data.frame(newdata))
  }
  
  
  n <- object@n
  
  # non conditional simulations
  sample_test <- sample.int(m,min(100,m))
  iter <- 1
  test <- FALSE
  yc <- matrix(0,nsim,m+n)
  #RMSE_save <- c() #-!!-#
  
  while((!test&iter<1000)){  
    a <- sample.int(m+n, size = nsim, replace = FALSE, prob = NULL)
    ya <- rnorm(nsim,0,1)
    if (requireNamespace("DiceKriging", quietly = TRUE)){
      Cab <- DiceKriging::covMat1Mat2(object@covariance, X1 = as.matrix(newdata[a,]), X2 = as.matrix(newdata), 
                                      nugget.flag = object@covariance@nugget.flag)/object@covariance@sd2
    }
    yc <- yc + Cab* (as.matrix(ya-diag(yc[,a]))%*%t(as.matrix(rep(1,m+n))))
    diag(yc[1:nsim,a]) <- ya
    
    if(iter%%20==0){
      Cemp <- var(yc[,sample_test]) 
      Ctheo <- DiceKriging::covMatrix(object@covariance, X = as.matrix(newdata[sample_test,]))[[1]]/
        object@covariance@sd2 
      RMSE_C <- sqrt(mean((Cemp-Ctheo)^2))  
      #RMSE_save <- c(RMSE_save,RMSE_C) #-!!-#
      if(RMSE_C < 0.1){test <- TRUE} 
    }
    
    iter <- iter+1
  }  
  rm(list=c("Cemp","Ctheo","ya" ,"Cab","a" ,"sample_test"))
  
  
  yc <- yc*sqrt(object@covariance@sd2*object@n/(object@n-object@p))
  
  if (cond){	
    if(type=="UK"){
      y.predict.trend <- F.newdata %*% object@trend.coef
      if (requireNamespace("DiceKriging", quietly = TRUE)){
        c.newdata <- DiceKriging::covMat1Mat2(object@covariance, X1 = object@X, X2 = newdata, 
                                              nugget.flag = object@covariance@nugget.flag)
      }
      
      rm(list=c("newdata"))	
      
      Tinv.c.newdata <- backsolve(t(object@T), c.newdata, upper.tri = FALSE)
      
      rm(list=c("c.newdata"))		
      
      y.predict.complement <- t(Tinv.c.newdata) %*% object@z
      y.predict <- y.predict.trend + y.predict.complement
      
      rm(list=c("y.predict.complement"))
      
      y.predict <- as.numeric(y.predict)
      
      ZD <- t(yc[,(m+1):(m+n)])
      if(!identical(object@noise.var,numeric(0))){
        ZD <- ZD + matrix(rnorm(nsim*n,0,1),nrow=n)*matrix(rep(sqrt(object@noise.var),nsim),ncol=nsim)
      }
      
      xtilde <- backsolve(t(object@T), ZD, upper.tri = FALSE)
      Betatilde <- solve(t(object@M)%*%object@M	)%*%t(object@M)%*%xtilde 
      
      rm(list=c("xtilde"))
      
      ztilde <-  backsolve(t(object@T), ZD -  object@F%*%Betatilde  , 
                           upper.tri = FALSE)
      
      rm(list=c("ZD"))
      
      ytilde.predict.complement <- t(Tinv.c.newdata) %*% ztilde
      
      rm(list=c("ztilde","Tinv.c.newdata"))
      
      ytilde.predict.trend <- F.newdata %*% Betatilde  
      
      rm(list=c("F.newdata","Betatilde"))
      
      ytilde.predict <- ytilde.predict.trend + ytilde.predict.complement
      
      rm(list=c("ytilde.predict.trend","ytilde.predict.complement","y.predict.trend"))
      
      yc <- t(matrix(rep(y.predict,nsim),ncol=nsim)) - t(ytilde.predict) +yc
      rm(list=c("y.predict","ytilde.predict"))
    }
    if(type=="SK"){
      y.predict.trend <- F.newdata %*% object@trend.coef
      if (requireNamespace("DiceKriging", quietly = TRUE)){
        c.newdata <- DiceKriging::covMat1Mat2(object@covariance, X1 = object@X, X2 = newdata, 
                                              nugget.flag = object@covariance@nugget.flag)
      }
      Tinv.c.newdata <- backsolve(t(object@T), c.newdata, upper.tri = FALSE)
      y.predict.complement <- t(Tinv.c.newdata) %*% object@z
      y.predict <- y.predict.trend + y.predict.complement
      y.predict <- as.numeric(y.predict)
      
      ZD <- t(yc[,(m+1):(m+n)])
      if(!identical(object@noise.var,numeric(0))){
        ZD <- ZD + matrix(rnorm(nsim*n,0,1),nrow=n)*matrix(rep(sqrt(object@noise.var),nsim),ncol=nsim)
      }
      
      ztilde <-  backsolve(t(object@T), ZD  ,  upper.tri = FALSE)
      ytilde.predict <- t(Tinv.c.newdata) %*% ztilde
      
      rm(list=c("ztilde","c.newdata","Tinv.c.newdata","newdata","y.predict.trend","F.newdata","ZD"))
      
      yc <- t(matrix(rep(y.predict,nsim),ncol=nsim)) - t(ytilde.predict) +yc
      rm(list=c("y.predict","ytilde.predict"))
    }
  }
  
  return(yc[,-c((m+1):(m+n))])  
  
}

##############################################################################

sobolpickfreeze <- function(y1,y2,nboot){
  output <- (apply(y1*y2,1,mean) - apply(y1,1,mean)*apply(y2,1,mean))/apply(y1,1,var)
  if(nboot > 1){
    n <- dim(y1)[2]
    for (i in 1:nboot){
      b <- sample(n,replace = TRUE)
      output <- rbind(output,(apply(y1[,b]*y2[,b],1,mean) - apply(y1[,b],1,mean)*apply(y2[,b],1,mean))/apply(y1[,b],1,var))
    }
  }
  return(output)
}

soboljansenpickfreeze <- function(y1,y2,nboot){
  varY <- apply(cbind(y1,y2),1,var)
  output <- (varY-apply((y1-y2)^2,1,mean)/2)/varY
  if(nboot > 1){
    n <- dim(y1)[2]
    for (i in 1:nboot){
      b <- sample(n,replace = TRUE)
      varY <- apply(cbind(y1[,b],y2[,b]),1,var)		
      output <- rbind(output,(varY-apply((y1[,b]-y2[,b])^2,1,mean)/2)/varY)
    }
  }
  return(output)
}

sobolEffpickfreeze <- function(y1,y2,nboot){
  output <- (apply(y1*y2,1,mean) - apply(cbind(y1,y2),1,mean)^2)/apply(cbind(y1,y2),1,var)
  if(nboot > 1){
    n <- dim(y1)[2]
    for (i in 1:nboot){
      b <- sample(n,replace = TRUE)
      output <- rbind(output,(apply(y1[,b]*y2[,b],1,mean) - apply(cbind(y1[,b],y2[,b]),1,mean)^2)/apply(cbind(y1[,b],y2[,b]),1,var))
    }
  }
  return(output)
}

sobolT2002pickfreeze <- function(y1,y2,nboot){
  output <- (1-(apply(y1*y2,1,mean) - apply(y1,1,mean)^2)/apply(y1,1,var))
  if(nboot > 1){
    n <- dim(y1)[2]
    for (i in 1:nboot){
      b <- sample(n,replace = TRUE)
      output <- rbind(output,(1-(apply(y1[,b]*y2[,b],1,mean) - apply(y1[,b],1,mean)^2)/apply(y1[,b],1,var)))
    }
  }
  return(output)
}

sobolTjansenpickfreeze <- function(y1,y2,nboot){
  output <- (apply((y1-y2)^2,1,mean)/2)/apply(y1,1,var)
  if(nboot > 1){
    n <- dim(y1)[2]
    for (i in 1:nboot){
      b <- sample(n,replace = TRUE)
      output <- rbind(output,(apply((y1[,b]-y2[,b])^2,1,mean)/2)/apply(y1[,b],1,var))
    }
  }
  return(output)
}

sobol2007pickfreeze <- function(y1,y2,y3,nboot){
  output <- apply(y3*(y2-y1),1,mean)/apply(y1,1,var)
  if(nboot > 1){
    n <- dim(y1)[2]
    for (i in 1:nboot){
      b <- sample(n,replace = TRUE)
      output <- rbind(output,apply(y3[,b]*(y2[,b]-y1[,b]),1,mean)/apply(y1[,b],1,var))
    }
  }
  return(output)
}

sobolT2007pickfreeze <- function(y1,y2,nboot){
  output <- apply(y1*(y1-y2),1,mean)/apply(y1,1,var)
  if(nboot > 1){
    n <- dim(y1)[2]
    for (i in 1:nboot){
      b <- sample(n,replace = TRUE)
      output <- rbind(output,apply(y1[,b]*(y1[,b]-y2[,b]),1,mean)/apply(y1[,b],1,var))
    }
  }
  return(output)
}
