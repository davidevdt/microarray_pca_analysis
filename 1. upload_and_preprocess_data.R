###########################################################################
###			UPLOAD AND PREPROCESS 2007 and 2008 TIV AND RMA DATA 		###
###########################################################################
#
# In this script, all the necessary datasets for the analysis will be  
# uploaded, prepocessed, and stored in the DATA folder. We will follow
# the following order: 
# 1. upload/process RMA 2008 dataset 
# 2. upload/process TIV 2008 outcome 
# 3. upload/process RMA 2007 dataset 
# 4. upload/process TIV 2007 outcome 
#
# Use setwd("<DIRECTORY>") to specify the working directory 
# --------------------------------------------------------------------------

######################################
#	0. Load Packages  and Gene Names #
######################################
library(dplyr)

# --------------------------------------------------------------------------
# Needed for loading RMA 2007 dataset and Gene Names:
# --------------------------------------------------------------------------
packages <- c("affy", "GEOquery")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(setdiff(packages, rownames(installed.packages())))  
}

packages <- c("Biobase", "annotate", "hgu133plus2.db")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
	source("http://bioconductor.org/biocLite.R")
	biocLite(setdiff(packages, rownames(installed.packages())))  
}

library(Biobase)
library(annotate)
library("hgu133plus2.db")
library(affy)
library(GEOquery)

# --------------------------------------------------------------------------
# Retrieve Gene Symbols (following https://github.com/katrijnvandeun/SPCovR/blob/master/R/Script_sgcca_spls.R )
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# Get and export official gene symbols: TAKES SOME TIME!
# --------------------------------------------------------------------------

x <- hgu133plus2SYMBOL 
symbolall<-Lkeys(x)
for (i in 1:length(symbolall)) {
	symbolall[i]<-get(symbolall[i],env=hgu133plus2SYMBOL)
}
write.table(symbolall, file="./DATA/symbolall.txt", col.names=FALSE, row.names=FALSE)

# --------------------------------------------------------------------------
# Get list of immune-system genes 
# --------------------------------------------------------------------------

imm_symb <- scan(file = "./DATA/list_immune_genes.txt", what = character())

# --------------------------------------------------------------------------

#######################
#	1. RMA 2008 data  #
#######################
# The dataset is already preprocessed and stored in the DATA folder (file TIVD3_rev.txt)
# --------------------------------------------------------------------------

#######################
#	2. TIV 2008 data  #
#######################
# We follow  https://github.com/katrijnvandeun/SPCovR/blob/master/MATLAB/ScriptHAI_TIVD28vsD0_2008.m (Matlab version)
# We will select the outcome at times D0
# --------------------------------------------------------------------------

y <- read.table("./DATA/TIVtiters.txt", header = FALSE, sep = "", na.strings = "NA")  

# --------------------------------------------------------------------------
# Subject ID 
# --------------------------------------------------------------------------

Subject=c(2,2,16,16,16,29,29,29,3,3,3,32,32,32,35,35,35,38,38,38,39,39,39,
			4,4,42,42,42,43,43,43,44,44,44,46,46,46,47,47,47,48,48,48,51,51,53,53,
			53,63,63,63,65,65,65,68,68,68,70,70,72,72,72,73,73,73,74,74,74,78,78,78,
			80,80,80,83,83,83,85,85,85)
			
# --------------------------------------------------------------------------
# Time ID 
# --------------------------------------------------------------------------

Time=c('D0','D7','D0','D3','D7','D0','D3','D7','D0','D3','D7','D0','D3','D7',
		'D0','D3','D7','D0','D3','D7','D0','D3','D7','D0','D3',
		'D0','D3','D7','D0','D3','D7','D0','D3','D7','D0','D3','D7',
		'D0','D3','D7','D0','D3','D7','D0','D3','D0','D3','D7','D0','D3',
		'D7','D0','D3','D7','D0','D3','D7','D0','D7','D0','D3','D7','D0','D3','D7',
		'D0','D3','D7','D0','D3','D7','D0','D3','D7','D0','D3', 'D7',
		'D0','D3','D7')			

# --------------------------------------------------------------------------		
# Create data frame and use only subjects with at least 2 time points
# --------------------------------------------------------------------------
dat <- data.frame(Subject=Subject, Time=Time)
n_times <- dat %>% group_by(Subject) %>% summarise(n_times = length(Time)) %>% filter(n_times > 2) 
dat <- dat %>% mutate(To_Retain = ifelse(Subject %in% n_times$Subject, 1, 0) )

# --------------------------------------------------------------------------
# D0
# --------------------------------------------------------------------------
matchD0 <- which(dat$Time == 'D0' & dat$To_Retain == 1)
Subj_d0 <- dat$Subject[matchD0]

# --------------------------------------------------------------------------
# Calculate antibody titers
# --------------------------------------------------------------------------
orig_titers <- read.table("./DATA/TIVtiters.txt", header = FALSE, sep = "", na.strings = "NA") 
sel_titers <- orig_titers[ , matchD0]
m1 <- sel_titers[2, ] / sel_titers[1,]
m2 <- sel_titers[4, ] / sel_titers[3,]
m3 <- sel_titers[6, ] / sel_titers[5,]
m <- rbind(m1, m2, m3)
M <- apply(m, 2, max)

# --------------------------------------------------------------------------
# Store TIV 2008 outcome 
# --------------------------------------------------------------------------
write.table(M, "./DATA/titer2008.txt", col.names = FALSE, row.names = FALSE)


#######################
#	3. RMA 2007 data  #
#######################
# We follow https://github.com/katrijnvandeun/SPCovR/blob/master/MATLAB/ScriptHAI_TIVD28vsD0_2007.m (Matlab version)
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# a. Get RMA 2007 data 
# --------------------------------------------------------------------------

accnr<-c("GSE29614") 			# Enter here the accession number: GSE29614 for TIV 2007

# --------------------------------------------------------------------------
# Load series and platform data from GEO
# --------------------------------------------------------------------------
getGEOSuppFiles( accnr )
setwd( paste("./", accnr,sep="") )
untar( paste(accnr,"_RAW.tar",sep="" ), exdir="rawdata" )
setwd( "./rawdata" )
cels <- list.celfiles()
data<-ReadAffy( filenames=cels )
# show(data): 27 samples = 9 subjects x 3 times (D0,D3,D7);

eset <- rma(data) 				# Create RMA pre-processed expression matrix

# --------------------------------------------------------------------------
# Store RMA 2007 dataset  
# --------------------------------------------------------------------------
setwd("../../")
write.exprs(eset, sep = "\t", file=("DATA/TIVRMA_2007.txt"), quote=FALSE, row.names = FALSE, col.names = FALSE)

# --------------------------------------------------------------------------
# b. Preprocess RMA 2007 data 
# --------------------------------------------------------------------------

# We follow  https://github.com/katrijnvandeun/SPCovR/blob/master/MATLAB/Script_TIVD3vsD0_2007.m (Matlab version)
tivrma <- read.table("./DATA/TIVRMA_2007.txt")

# --------------------------------------------------------------------------
# Subject ID 
# --------------------------------------------------------------------------
Subject = c(12, 12, 12, 16, 16, 16, 18, 18, 18, 21, 21, 21, 23, 23, 23, 29, 29, 29, 32, 32, 32, 35, 35, 35, 39, 39, 39)

# --------------------------------------------------------------------------
# Time ID 
# --------------------------------------------------------------------------
Time = c('D0','D3','D7','D0','D3','D7',
		 'D0','D3','D7','D0','D3','D7',
		 'D0','D3','D7','D0','D3','D7',
		 'D0','D3','D7','D0','D3','D7',
		 'D0','D3','D7')

dat <- data.frame(Subject = Subject, Time = Time)		 
		 
# --------------------------------------------------------------------------
# D0	
# --------------------------------------------------------------------------
foundD0 <- which(dat$Time == 'D0')
Subj_d0 <- dat$Subject[foundD0]

# --------------------------------------------------------------------------
# D3	
# --------------------------------------------------------------------------
foundD3 <- which(dat$Time == 'D3')
Subj_d3 <- dat$Subject[foundD3]


TIVBLOCKD0 <- tivrma[ ,foundD0]
TIVBLOCKD3 <- tivrma[ ,foundD3]
TIVBLOCKD3 <- TIVBLOCKD3 - TIVBLOCKD0 

# --------------------------------------------------------------------------
# Scale (transpose of) data 
# --------------------------------------------------------------------------
TIVBLOCKD3_std <- scale( t(TIVBLOCKD3) )

# --------------------------------------------------------------------------
# Select immune genes 
# --------------------------------------------------------------------------

imm_ind <- which(symbolall %in% imm_symb)
# length(imm_ind)
TIVBLOCKD3_std <- TIVBLOCKD3_std[,imm_ind]	
colnames(TIVBLOCKD3_std) <- make.names(colnames(TIVBLOCKD3_std), unique=T) 		# We differentiate genes with the same name
rownames(TIVBLOCKD3_std) <- 1:nrow(TIVBLOCKD3_std)

# --------------------------------------------------------------------------
# Store RMA 2007 (preprocessed) dataset 
# --------------------------------------------------------------------------
write.table(TIVBLOCKD3_std, file = "./DATA/TIVD3_2007_rev.txt")
		 
#######################
#	4. TIV 2007 data  #
#######################
# We follow https://github.com/katrijnvandeun/SPCovR/blob/master/MATLAB/ScriptHAI_TIVD28vsD0_2007.m (Matlab version)
# We will select the outcome at times D0
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# Subject ID 
# --------------------------------------------------------------------------
Subject = c(12, 12, 12, 16, 16, 16, 18, 18, 18, 21, 21, 21, 23, 23, 23, 29, 29, 29, 32, 32, 32, 35, 35, 35, 39, 39, 39)

# --------------------------------------------------------------------------
# Time ID 
# --------------------------------------------------------------------------
Time = c('D0','D3','D7','D0','D3','D7',
		 'D0','D3','D7','D0','D3','D7',
		 'D0','D3','D7','D0','D3','D7',
		 'D0','D3','D7','D0','D3','D7',
		 'D0','D3','D7')

# --------------------------------------------------------------------------		
# Create data frame and use only subjects with at least 2 time points
# --------------------------------------------------------------------------
dat <- data.frame(Subject=Subject, Time=Time)
n_times <- dat %>% group_by(Subject) %>% summarise(n_times = length(Time)) %>% filter(n_times > 2) 
dat <- dat %>% mutate(To_Retain = ifelse(Subject %in% n_times$Subject, 1, 0) )

# --------------------------------------------------------------------------
# D0
# --------------------------------------------------------------------------
matchD0 <- which(dat$Time == 'D0' & dat$To_Retain == 1)
Subj_d0 <- dat$Subject[matchD0]

# --------------------------------------------------------------------------
# Calculate antibody titers
# --------------------------------------------------------------------------
orig_titers <- read.table("./DATA/TIVtiters_2007.txt", header = FALSE, sep = "", na.strings = "NA") 
sel_titers <- orig_titers[ , matchD0]
m1 <- sel_titers[2, ] / sel_titers[1,]
m2 <- sel_titers[4, ] / sel_titers[3,]
m3 <- sel_titers[6, ] / sel_titers[5,]
m <- rbind(m1, m2, m3)
M <- apply(m, 2, max)

# --------------------------------------------------------------------------
# Store TIV 2007 outcome 
# --------------------------------------------------------------------------
write.table(M, "./DATA/titer2007.txt", col.names = FALSE, row.names = FALSE)

