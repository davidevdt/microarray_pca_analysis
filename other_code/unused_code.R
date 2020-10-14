################################# 
### spca: plot weights matrix ###
#################################
# --------------------------------------------------------------------------
# Rearrange weights to favour block structure visualization: 
# --------------------------------------------------------------------------
Wmat_transformed <- Wmat
Wmat_transformed[Wmat==0] <- -Inf
Wmat_transformed <- data.frame(Wmat_transformed, colnames(X), srt = 1:nrow(Wmat))
names(Wmat_transformed) <- c("Component.1", "Component.2", "rownames.Wmat.", "srt")
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

# --------------------------------------------------------------------------
### Plot results 
# --------------------------------------------------------------------------
base_size <- 9
brks <- seq(1,ncol(X)*num_comp, by = 10)
WmatVec <- c(as.matrix(new_Wmat))
VarNames <- factor(rep(new_var_names, num_comp))
lbls <- as.character(VarNames) 
lbls[-brks] <- ""
Comps <- rep(paste0("Component ", 1:num_comp), each = ncol(X))
WMatDataFrame <- data.frame(Weights = WmatVec, Variables = VarNames, Component = Comps, Sort = all_sorted)

p3 <- ggplot(WMatDataFrame, aes(x = Component, y = reorder(Variables, sort(Sort, decreasing = TRUE)))) + 
		geom_tile(aes(fill = Weights)) + 
 		scale_fill_gradient2(low = "red", mid = "black", high = "green", limits = c(-0.02,0.02), oob = scales::squish) +
 		labs(y = "Genes") + 
 		scale_y_discrete(labels=lbls)
		
p3 + theme_grey(base_size = base_size) + labs(x = "", y = "") +
 				scale_x_discrete(expand = c(0, 0)) 
				
p3 + theme(axis.text.x = element_text(size=10, color = "black"),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(),
            legend.title = element_text(colour="black", size=10),
            legend.text = element_text(colour="black", size=10))
	

############################# 
### vbpca: unsorted plots ###
#############################

# --------------------------------------------------------------------------
### Plot results (unsorted)
# --------------------------------------------------------------------------

plotUnsorted <- TRUE
plotProb <- FALSE 
bounds <- c(-0.1, 0.1)

if (plotUnsorted){ 
	if (plotProb) {
	
		# Probability Matrix 
		ProbMatVec <- c(bpca_mod[[9]])
		VarNames <- factor(rep(colnames(X), num_comp))
		lbls <- as.character(VarNames) 
		lbls[-brks] <- ""
		Comps <- rep(paste0("Component ", 1:num_comp), each = ncol(X))
		PrMatDataFrame <- data.frame(Probability = ProbMatVec, Variables = VarNames, Component = Comps)
		p4 <- ggplot(PrMatDataFrame, aes(x = Component, y = VarNames)) + 
				geom_tile(aes(fill = Probability)) + 
				scale_fill_gradient(low = "white", high = "blue") + 
				labs(y = "Genes") + 
				scale_y_discrete(labels=lbls)

		p4 + theme_grey(base_size = base_size) + labs(x = "", y = "") +
						scale_x_discrete(expand = c(0, 0)) 
					
		p4 + theme(axis.text.x = element_text(size=10, color = "black"),
					 axis.text.y = element_text(size = 10, color = "black"),
					 axis.title.x = element_blank(),
					 axis.title.y = element_blank(),
					 axis.ticks.x=element_blank(),
					 axis.ticks.y=element_blank(),
					 legend.title = element_text(colour="black", size=10),
					 legend.text = element_text(colour="black", size=10))
		} else {

			#  Weights Matrix 
			WmatVec <- c(Wmat)
			WMatDataFrame <- data.frame(Weights = WmatVec, Variables = VarNames, Component = Comps)

			p5 <- ggplot(WMatDataFrame, aes(x = Component, y = VarNames)) + 
					geom_tile(aes(fill = Weights)) + 
					scale_fill_gradient2(low = "red", mid = "black", high = "green", limits = bounds,  oob = scales::squish) +
					labs(y = "Genes") + 
					scale_y_discrete(labels=lbls)

			p5 + theme_grey(base_size = base_size) + labs(x = "", y = "") +
							scale_x_discrete(expand = c(0, 0)) 
			 p5 + theme(axis.text.x = element_text(size=10, color = "black"),
						 axis.text.y = element_text(size = 10, color = "black"),
						 axis.title.x = element_blank(),
						 axis.title.y = element_blank(),
						 axis.ticks.x=element_blank(),
						 axis.ticks.y=element_blank(),
						 legend.title = element_text(colour="black", size=10),
						 legend.text = element_text(colour="black", size=10))
		}
}
	