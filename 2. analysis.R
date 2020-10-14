####################################################
###			ANALYZE 2008 and 2007 DATA  		###
###################################################
# 
# In this script, we will perform model selection for the PCA models (elasticnet Sparse PCA and Variational Bayes PCA)
# Model selection will be performed on the 2008 RMA dataset 
# Subsequently, we will perform PCA regression using the  2008 and 2007 RMA datasets to define the principal components, 
# and the corresponding 2008 and 2007 TIV outcomes as outcomes. 
# 
# We will proceed as follows: 
# 1. PCA model selection (number of components - spca model selection - vbpca model selection) on 2008 RMA dataset 
# 2. PCA regression using 2008 data 
# 3. PCA regression using 2007 data 
# 
# Use setwd("<DIRECTORY>") to specify the working directory 
# --------------------------------------------------------------------------

#################################################
#	0. Load Packages and retrieve Gene Symbols  #
#################################################
library(dplyr)
# devtools::install_github("davidevdt/bayespca")
library(bayespca)
library(ggplot2)
source("./eigencv.R")

# --------------------------------------------------------------------------
# Load data
# --------------------------------------------------------------------------

X <- read.table("./DATA/TIVD3_rev.txt", header = FALSE, sep = "", na.strings = "NA")
X <- as.matrix(X)
X <- t(X)
X <- X[ ,1:54675]
rownames(X) <- 1:nrow(X)
# Data is already scaled : round( apply(X, 2, mean), 1 ) and round( apply( X, 2, sd), 1 )


# --------------------------------------------------------------------------
# Upload gene names and immune genes 
# --------------------------------------------------------------------------
symbolall <- scan(file = "./DATA/symbolall.txt", what = character())
imm_symb <- scan(file = "./DATA/list_immune_genes.txt", what = character())

# --------------------------------------------------------------------------
# Select immune genes
# --------------------------------------------------------------------------
imm_ind <- which(symbolall %in% imm_symb)
# length(imm_ind)
X <- X[,imm_ind]										# X is going to be the working dataset (2008 RMA data)
colnames(X) <- symbolall[imm_ind]
copy_names <- colnames(X)								# We store the name of the original genes
colnames(X) <- make.names(colnames(X), unique=T) 		# We differentiate genes with the same name


#############################
#	1. PCA model selection  #
#############################
# --------------------------------------------------------------------------
# a. Number of Principal Components (with SVD)
# --------------------------------------------------------------------------
svd_mod <- svd(X)	
s_vals <- svd_mod$d
barplot((s_vals^2/sum(s_vals^2))[1:10], names.arg = 1:10, col = "blue", lwd = 2, main = "", 
				xlab = "N. Comp", ylab = "Variance Accounted For")	

cum_perc_variance <- cumsum(svd_mod$d^2) / sum(svd_mod$d^2)
cum_perc_variance
plot(cum_perc_variance, type="b", xlab="N. Comp", ylab="Cumulative % variance", col="blue")

# --------------------------------------------------------------------------
# We will evaluate the first 2 components (~27.43% of total variance). 
# --------------------------------------------------------------------------
num_comp <- 2

# --------------------------------------------------------------------------
# b. Regularized PCA - spca 
# --------------------------------------------------------------------------
#
# Searching for the max. penalty values which leads to only 1 non-zero weight: 

tol <- 1e-04
grid_length <- 20
max_iter <- 1e+05

min_alpha <- 1e-5
max_penalty <- 700
grid_A <- seq(min_alpha, max_penalty, length.out = grid_length )
i <- 1 

while( TRUE ){
	cat("Penalty = ", grid_A[i], " , ", i, " of ", length(grid_A), "\n")
	tmp_mod <- elasticnet::spca(x = X, K = num_comp, para = rep(grid_A[i], num_comp), sparse = "penalty", max.iter = max_iter, 
								   eps.conv = tol, trace = FALSE, 
								   type = "predictor", use.corr = FALSE)
							   
	if( sum(tmp_mod$loadings) == 1 ) {
		cat("Found: max penalty = ", grid_A[i], "\n" )
		max_alpha <- grid_A[i]
		break 
	}
			
	i = i + 1

}

# --------------------------------------------------------------------------
# Now let's tune the SPCA model with lambda large (lambda --> Infinity) 	
# --------------------------------------------------------------------------

lambda <- 1e+6
grid_length <- 20
min_alpha <- 1e-5
grid_B <- (seq((min_alpha), (max_alpha), length.out = grid_length) )		# Lasso penalty grid
numFolds <- 5																

all_Mse <- rep(0, grid_length)

# Cross-Validation with S.E. rule: 
mse <- rep(0, length(grid_B))
minMse <- Inf
sdMse <- 0 
selPar <- NULL			
			
for( i in 1:length(grid_B)  ){
	cat("i = ", i, " of ", length(grid_B), " (", grid_B[i] ,")", "\n")				
	set.seed(71)
	eCV <- EigenCVel_net(X, alpha = grid_B[i], 
							beta = lambda, nFolds = numFolds, 
							D = num_comp, maxIt = max_iter, tol = tol, 
							scaleDat = FALSE, verbose = FALSE)
	mse[i] <- eCV$MSE
	all_Mse[i] <- mse[i]
				
	if(mse[i] < minMse){
		minMse <- mse[i]
		sdMse <- eCV$sdMSE
	}	
								
}
			
indx <- which(mse <= minMse + sdMse )
selPar <- max(indx)
cat("Selected lasso parameter: ", grid_B[selPar], ", index=", selPar, "\n" )
alphaSel <- grid_B[selPar]	

# --------------------------------------------------------------------------
# Plot CV results 
# --------------------------------------------------------------------------

plot(grid_B, all_Mse, main="CV Results", type="b", col=ifelse(all_Mse <= (min(all_Mse) + sdMse),"white","blue"), 
											xlab="Penalty", ylab="MSE", xaxt="n", lwd=2, 
											ylim=c(max(0,(min(all_Mse)) - sdMse), max(all_Mse)))
											
lines(grid_B, rep((min(all_Mse) + sdMse), length(grid_B)), col="blue", lty=2, lwd=2)

lines(grid_B, rep((max(0,(min(all_Mse)) - sdMse)), length(grid_B)), col="blue", lty=2, lwd=2)

points(grid_B[which(all_Mse <= (min(all_Mse) + sdMse))], 
			all_Mse[which(all_Mse <=  (min(all_Mse) + sdMse))], col="red", lwd=2, pch="x", cex=2)
			
points(grid_B[max(which(all_Mse <= (min(all_Mse) + sdMse)))], 
		all_Mse[max(which(all_Mse <=  (min(all_Mse) + sdMse)))], col="green", lwd=2, pch="x", cex=2)
		
axis(1, at = grid_B, las=1)

# --------------------------------------------------------------------------
# Model Estimation with Selected Hyperparameters 
# --------------------------------------------------------------------------
mod <- elasticnet::spca(x = X, K = num_comp, para = rep(alphaSel, num_comp), 
						lambda = lambda, sparse = "penalty", 
						max.iter = max_iter, eps.conv = tol, 
						trace = FALSE, type = "predictor", use.corr = FALSE)	
						
# % of non-zero's in PC-1						
j <- 1 
sum(mod$loadings[,j]!=0) / nrow(mod$loadings)			

# % of non-zero's in PC-2
j <- 2
sum(mod$loadings[,j]!=0) / nrow(mod$loadings)			

# Overall % of non-zero's  		   
sum(mod$loadings!=0) / (nrow(mod$loadings)	* num_comp)			

# --------------------------------------------------------------------------
# Reconstruction error 
# --------------------------------------------------------------------------
Wmat <- mod$loadings
XTX <- t(X) %*% X
sVd <- svd(XTX %*% mod$loadings)
P <- sVd$u %*% sVd$v		
rec_err_spca <- sum( (X - ( X %*% mod$loadings %*% t(P) ) )^2)			# spca --> Reconstruction Error 

# --------------------------------------------------------------------------
# c. Regularized PCA - vbpca 
# --------------------------------------------------------------------------
#
# Instead of heuristically picking the hyperparameters alpha and beta of the InverseGamma priors, 
# we will use a data-driven procedure: we will try different values of the hyperparameters
# and will pick the ones which lead to largest ELBO
# (Notice that this does not require splitting the dataset into folds like in Cross-Validations)
# 
# --------------------------------------------------------------------------

alphatau_grid <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20)				# Grid for alpha 
betatau_grid <- c(1, 2, 3, 5, 6, 7, 8, 9, 10)							# Grid for beta 

vbpcaPars <- expand.grid(alpha = alphatau_grid, beta = betatau_grid)
vbpcaPars$recerr <- rep(0, nrow(vbpcaPars))
vbpcaPars$nonZeroLoadings <- rep(0, nrow(vbpcaPars))
vbpcaPars$elbo <- rep(0, nrow(vbpcaPars))

max_iter <- 1e+05
tol <- 1e-04															# Convergence tolerance 
useSimplifiedWMat <- TRUE												# Use the weight matrix with un-selected weights set to 0 
v0_par <- 1e-04															# v parameter for spike variance 

for( i in 1:nrow(vbpcaPars)  ){
	cat("i = ", i, " of ", nrow(vbpcaPars), "\n")
	
 	ctrl <- vbpca_control(center = FALSE, scalecorrection = -1, alphatau = vbpcaPars$alpha[i], 
							betatau = vbpcaPars$beta[i], gammatau = -1, hypertype = "common", 
							beta1pi = 1, beta2pi = 1, v0 = v0_par, plot.lowerbound = FALSE)
							
 	bpca_mod <- vbpca(X = X, D = num_comp, maxIter = max_iter, verbose = FALSE, 
						priorvar = 'invgamma', updatetau = FALSE, global.var = FALSE,
						SVS = TRUE, priorInclusion = 0.5, tolerance = tol,
						control = ctrl, suppressWarnings = TRUE)
	
 	if ( !useSimplifiedWMat ) {
 		recErr <- sum( (X - ( X %*% bpca_mod[[1]] %*% t(bpca_mod[[2]]) ) )^2)								# Reconstruction Error 	
 		vpcaPars$recerr[i] <- recErr																		
		vbpcaPars$nonZeroLoadings[i] <-  sum(bpca_mod[[9]]>=0.5) / (nrow(bpca_mod[[9]]) * num_comp)			# % non-zero weights 
		vbpcaPars$elbo[i] <- bpca_mod$elbo 																	# Model ELBO 
 	} else {
 		ind <- (bpca_mod[[9]] < 0.5)			
 		Wmat <- bpca_mod[[1]]
 		Wmat[ind == TRUE] <- 0 
 		recErr <- sum( (X - ( X %*% Wmat %*% t(bpca_mod[[2]]) ) )^2)										# Reconstruction Error 	
 		vbpcaPars$recerr[i] <- recErr
		vbpcaPars$nonZeroLoadings[i] <-  sum(bpca_mod[[9]]>=0.5) / (nrow(bpca_mod[[9]]) * num_comp)			# % non-zero weights 
		vbpcaPars$elbo[i] <- bpca_mod$elbo 																	# Model ELBO 
 	}
}

# --------------------------------------------------------------------------
# Filter models with % non-zero weights == 0 or 1 
# --------------------------------------------------------------------------

vbpcaPars2 <- vbpcaPars[vbpcaPars$nonZeroLoadings != 1 & vbpcaPars$nonZeroLoadings != 0, ]

# --------------------------------------------------------------------------
# Select model with largest ELBO 
# --------------------------------------------------------------------------
sel_pars <- vbpcaPars2[which.max(vbpcaPars2$elbo),]
sel_alpha <- sel_pars$alpha 
sel_beta <- sel_pars$beta 
cat("Selected InverseGamma hyperparameters: alpha = ", sel_alpha,  " beta =  ", sel_beta, "\n" )

# --------------------------------------------------------------------------
# Estimate the chosen model
# --------------------------------------------------------------------------  
ctrl <- vbpca_control(center = FALSE, scalecorrection = -1, alphatau = sel_alpha, betatau = sel_beta, 
						beta1pi = 1, beta2pi = 1, v0 = v0_par)

bpca_mod <- vbpca(X = X, D = num_comp, maxIter = max_iter, verbose = TRUE, 
					priorvar = 'invgamma', updatetau = FALSE, global.var = FALSE,
					SVS = TRUE, priorInclusion = 0.5, tolerance = tol,
					control = ctrl)

# % of non-zero's in PC-1						
j <- 1 
sum(bpca_mod[[9]][,j] >= 0.5) / nrow(bpca_mod[[9]])	

# % of non-zero's in PC-2
j <- 2
sum(bpca_mod[[9]][,j] >= 0.5) / nrow(bpca_mod[[9]])		

# Overall % of non-zero's  		   
sum(bpca_mod[[9]] >= 0.5) / (nrow(bpca_mod[[9]]) * num_comp)

# --------------------------------------------------------------------------
# Reconstruction error 
# --------------------------------------------------------------------------
ind <- (bpca_mod[[9]] <= 0.5)			
Wmat <- bpca_mod[[1]]
Wmat[ind == TRUE] <- 0 	
rec_err_vbpca <- sum( (X - ( X %*% Wmat %*% t(bpca_mod[[2]]) ) )^2)			# vbpca --> Reconstruction Error 


# --------------------------------------------------------------------------
# bpca vs spca Reconstruction errors 
# --------------------------------------------------------------------------
c("Rec_Error_vbpca" = rec_err_vbpca, "Rec_Error_spca" = rec_err_spca)


# --------------------------------------------------------------------------
# Plot weights and probabilities heatmaps (vbpca)
# --------------------------------------------------------------------------
#
# Rearrange weights to favour block structure visualization: 

Wmat_transformed <- Wmat
Wmat_transformed[Wmat==0] <- -Inf
Wmat_transformed <- data.frame(Wmat_transformed, rownames(Wmat), srt = 1:nrow(Wmat))

ProbMat <- data.frame(bpca_mod[[9]]) %>% arrange(desc(Wmat_transformed$Component.1), desc(Wmat_transformed$Component.2))

Wmat_transformed <- Wmat_transformed %>% arrange(desc(Component.1), desc(Component.2))

new_Wmat <- Wmat_transformed[,1:2]
new_Wmat[new_Wmat == -Inf] <- 0

new_var_names <- Wmat_transformed$rownames.Wmat.
new_var_names <- new_var_names[length(new_var_names): 1]

new_sort <- Wmat_transformed$srt 

sorted_new_sort <- new_sort[length(new_sort):1]

lst <- vector("list", num_comp)
for( nc in 1:num_comp  ){
 	if(nc == 1){
 		lst[[nc]] <- sorted_new_sort 
 	}else{
 		lst[[nc]] <- sorted_new_sort + max(lst[[nc-1]])
 	}
 }
all_sorted <- do.call("c", lst)

### Plotting options: 

base_size <- 9
brks <- seq(1,ncol(X)*num_comp, by = 50)

# Probability Matrix: 

ProbMatVec <- c(as.matrix(ProbMat))

VarNames <- factor(rep(new_var_names, num_comp))
lbls <- as.character(VarNames) 
lbls[-brks] <- ""

Comps <- rep(paste0("Component ", 1:num_comp), each = ncol(X))

PrMatDataFrame <- data.frame(Probability = ProbMatVec, Variables = VarNames, Component = Comps, Sort = all_sorted)

p1 <- ggplot(PrMatDataFrame, aes(x = Component, y = reorder(Variables, sort(Sort, decreasing = TRUE)) )) + 
 		geom_tile(aes(fill = Probability)) + 
 		scale_fill_gradient(low = "gray85", high = "blue") + 
 		labs(y = "Genes") + 
 		scale_y_discrete(labels=lbls)

p1 + theme_grey(base_size = base_size) + labs(x = "", y = "") +
				scale_x_discrete(expand = c(0, 0)) 
			
p1 + theme(axis.text.x = element_text(size=10, color = "black"),
           axis.text.y = element_text(size = 10, color = "black"),
           axis.title.x = element_blank(),
           axis.title.y = element_blank(),
           axis.ticks.x=element_blank(),
           axis.ticks.y=element_blank(),
           legend.title = element_text(colour="black", size=10),
           legend.text = element_text(colour="black", size=10))

# Weights Matrix : 

bounds <- c(-0.1, 0.1)
WmatVec <- c(as.matrix(new_Wmat))
WMatDataFrame <- data.frame(Weights = WmatVec, Variables = VarNames, Component = Comps, Sort = all_sorted)

p2 <- ggplot(WMatDataFrame, aes(x = Component, y = reorder(Variables, sort(Sort, decreasing = TRUE)) )) + 
				geom_tile(aes(fill = Weights)) + 
				scale_fill_gradient2(low = "red", mid = "black", high = "green", limits = bounds, oob = scales::squish) +
				labs(y = "Genes") + 
				scale_y_discrete(labels=lbls)
		
p2 + theme_grey(base_size = base_size) + labs(x = "", y = "") +
				scale_x_discrete(expand = c(0, 0)) 
				
 p2 + theme(axis.text.x = element_text(size=10, color = "black"),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(),
            legend.title = element_text(colour="black", size=10),
            legend.text = element_text(colour="black", size=10))
			
			
# --------------------------------------------------------------------------
# Export lists of genes identified by Baysian PCA 
# --------------------------------------------------------------------------		

genes_PC1 <- copy_names[(bpca_mod[[9]][,1] >= 0.5)] 
genes_PC2 <- copy_names[(bpca_mod[[9]][,2] >= 0.5)] 

write.table(genes_PC1, file="./DATA/genes_PC1.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(genes_PC2, file="./DATA/genes_PC2.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

# These lists will be exported in the page https://david.ncifcrf.gov/summary.jsp for functional annotation analysis. 
# (Select the Identifier: OFFICIAL_GENE_SYMBOL)

###################################
#	2. PCA regression - TIV 2008  #
###################################

# --------------------------------------------------------------------------
# Load Y (2008), transform it to log-scale, and center it 
# --------------------------------------------------------------------------
titers <-read.table("./DATA/titer2008.txt", header = FALSE, sep = "", na.strings = "NA", col.names = "")[,1]
y <- log2(titers)
y <- y - mean(y)

# --------------------------------------------------------------------------
# Create predictor matrix for Bayes PCA
# --------------------------------------------------------------------------
mat <- bpca_mod[[1]]
ind_zeros <- which(bpca_mod[[9]] <= 0.5)
Wmat[ind_zeros] <- 0 

Z_bpca <- data.frame(X %*% Wmat)
names(Z_bpca) = c("Component_1", "Component_2")

# --------------------------------------------------------------------------
# Linear Model with Z_bpca 
# --------------------------------------------------------------------------
mod_bpca <- lm(Y ~ ., data.frame(Y = y, Component_1 = Z_bpca$Component_1, Component_2 = Z_bpca$Component_2))
summary(mod_bpca)

# --------------------------------------------------------------------------
# Create predictor matrix for Sparse PCA
# --------------------------------------------------------------------------
Wmat_spca <- mod$loadings
Z_spca <- data.frame(X %*% Wmat_spca)
names(Z_spca) = c("Component_1", "Component_2")

# --------------------------------------------------------------------------
# Linear Model with Z_spca 
# --------------------------------------------------------------------------
mod_spca <- lm(Y ~ ., data.frame(Y = y, Component_1 = Z_spca$Component_1, Component_2 = Z_spca$Component_2))
summary(mod_spca)

# --------------------------------------------------------------------------
# R^2 for the two models 
# --------------------------------------------------------------------------
c("R_squared_bpca_2008" = summary(mod_bpca)$r.squared, "R_squared_spca_2008" = summary(mod_spca)$r.squared)

###################################
#	3. PCA regression - TIV 2007  #
###################################

# --------------------------------------------------------------------------
# Load X (2007) 
# --------------------------------------------------------------------------
X_2007 <- read.table("./DATA/TIVD3_2007_rev.txt", header = TRUE, sep = "", na.strings = "NA")
colnames(X_2007) <- symbolall[imm_ind]

# --------------------------------------------------------------------------
# Load Y (2007), transform it to log-scale, and center it 
# --------------------------------------------------------------------------
titers <-read.table("./DATA/titer2007.txt", header = FALSE, sep = "", na.strings = "NA", col.names = "")[,1]
y <- log2(titers)
y <- y - mean(y)

# --------------------------------------------------------------------------
# Create predictor matrix for Bayes PCA
# --------------------------------------------------------------------------
Wmat <- bpca_mod[[1]]
ind_zeros <- which(bpca_mod[[9]] <= 0.5)
Wmat[ind_zeros] <- 0 

Z_bpca <- data.frame(as.matrix(X_2007) %*% Wmat)
names(Z_bpca) = c("Component_1", "Component_2")

# --------------------------------------------------------------------------
# Linear Model with Z_bpca 
# --------------------------------------------------------------------------
mod_bpca_2007 <- lm(Y ~ ., data.frame(Y = y, Component_1 = Z_bpca$Component_1, Component_2 = Z_bpca$Component_2))
summary(mod_bpca_2007)

# --------------------------------------------------------------------------
# Create predictor matrix for Sparse PCA
# --------------------------------------------------------------------------
Wmat_spca <- mod$loadings

Z_spca <- data.frame(as.matrix(X_2007)  %*% Wmat_spca)
names(Z_spca) = c("Component_1", "Component_2")

# --------------------------------------------------------------------------
# Linear Model with Z_spca 
# --------------------------------------------------------------------------
mod_spca_2007 <- lm(Y ~ ., data.frame(Y = y, Component_1 = Z_spca$Component_1, Component_2 = Z_spca$Component_2))
summary(mod_spca_2007)

# --------------------------------------------------------------------------
# R^2 for the two models 
# --------------------------------------------------------------------------
c("R_squared_bpca_2007" = summary(mod_bpca_2007)$r.squared, "R_squared_spca_2007" = summary(mod_spca_2007)$r.squared)