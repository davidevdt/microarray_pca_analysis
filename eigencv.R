# EigenVector cross-validation, Bro et al. (2008)
EigenCV <- function(X, alpha, beta, nFolds, D, maxIt, tol, scaleDat, verbose){
	
	I <- nrow(X)
	J <- ncol(X)
	
	folds <- rep_len(1:nFolds, I)
	cvErr <- matrix(NA, I, J)
	errFold <- rep(NA, nFolds)
	
	
	for( f in 1:nFolds ){
		
		# Split of the data 
		currFold <- which(folds == f)
		Xtrain <- X[-currFold, ]
		Xtest <- X[currFold,]
		
		if(scaleDat){
			sc <- scaleX(Xtrain, center = TRUE, scalingFactor = 0)
			Xtrain <- sc
			attr(Xtrain, "center") <- NULL
			attr(Xtrain, "scaling") <- NULL
			
			Xtest <- t( t(Xtest) - attr(sc, "center"))
			Xtest <- t( t(Xtest) - attr(sc, "scaling"))
			rm(sc)
		}
	
		# Model estimation 
		mod <- sparsepca::spca(Xtrain, k = D, alpha = alpha, beta = beta, 
							   center = FALSE, scale = FALSE, max_iter = maxIt, 
							   tol = tol, verbose = FALSE)
	
		# EigenVector CV
		predicted <- matrix(NA, nrow(Xtest), J)
		if(verbose) cat("Fold: ", f, "\n")
		for(j in 1:J){
			scoresNoJ <- Xtest[,-j] %*% mod$loadings[-j, ]
			predicted[,j] <- scoresNoJ %*% mod$transform[j, ]	
		}
		cvErr[currFold, ] <- (Xtest - predicted)^2 
		errFold[f] <- mean(cvErr[currFold,])
		
	}
	
	return(list(
		cvError = cvErr, 
		MSE = mean(errFold), 
		MSEFold = errFold, 
		sdMSE = sd(errFold) / sqrt(nFolds)
	))	
}












# EigenCV for elasticnet package 
EigenCVel_net <- function(X, alpha, beta, nFolds, D, maxIt, tol, scaleDat, verbose){
	
	I <- nrow(X)
	J <- ncol(X)
	
	folds <- rep_len(1:nFolds, I)
	cvErr <- matrix(NA, I, J)
	errFold <- rep(NA, nFolds)
	
	
	for( f in 1:nFolds ){
		
		# Split of the data 
		currFold <- which(folds == f)
		Xtrain <- X[-currFold, ]
		Xtest <- X[currFold,]
		
		if(scaleDat){
			sc <- scaleX(Xtrain, center = TRUE, scalingFactor = 0)
			Xtrain <- sc
			attr(Xtrain, "center") <- NULL
			attr(Xtrain, "scaling") <- NULL
			
			Xtest <- t( t(Xtest) - attr(sc, "center"))
			Xtest <- t( t(Xtest) - attr(sc, "scaling"))
			rm(sc)
		}
	
		# Model estimation 
		mod <- elasticnet::spca(Xtrain, K = D, para = rep(alpha, D), lambda = beta,
								sparse = "penalty", max.iter = maxIt, 
								eps.conv = tol, trace = FALSE, 
								type = "predictor", use.corr = FALSE)
							   
							   
		XTX <- t(Xtrain) %*% Xtrain
		sVd <- svd(XTX %*% mod$loadings)
		P <- sVd$u %*% sVd$v
	
		# EigenVector CV
		predicted <- matrix(NA, nrow(Xtest), J)
		if(verbose) cat("Fold: ", f, "\n")
		for(j in 1:J){
			scoresNoJ <- Xtest[,-j] %*% mod$loadings[-j, ]
			predicted[,j] <- scoresNoJ %*% P[j, ]	
		}
		cvErr[currFold, ] <- (Xtest - predicted)^2 
		errFold[f] <- mean(cvErr[currFold,])
		
	}
	
	return(list(
		cvError = cvErr, 
		MSE = mean(errFold), 
		MSEFold = errFold, 
		sdMSE = sd(errFold) / sqrt(nFolds)
	))	
}

