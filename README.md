# Microarray PCA Analysis 

Prediction via PCA regression of the efficacy of vaccine against influenza. The micro-array gene expression data (2007 and 2008) are publicly available on the NCBI Gene Expression Omnibus
(https://www.ncbi.nlm.nih.gov/geo/ - accession numbers GSE29614 and GSE29617). The goal is to predict efficacy of the vaccine with gene expression data obtained after vaccination. <br/>
To achieve that, we perform PCA regression where Sparse PCA (```elasticnet``` package) and Bayesian PCA ([```bayespca```](https://github.com/davidevdt/bayespca) package ) are compared. Data retrieval and preprocessing follows from [2]. 

## Running the Analysis 
The analysis was carried out with R 3.5.0. Some of the datasets are alread included in the folder DATA. To replicate the analysis, perform the following operations: 

1. Clone the repository in a directory <DIRECTORY>
2. Open R and set the working directory as <DIRECTORY>
3. Launch the script "1. upload_and_prepocess_data.R" to upload the remaining datasets and obtain the gene symbols
4. Launch the script "2. analysis.R" to perform model selection and observe results of the analysis 

Points (3) and (4) (especially the model selection step of the Sparse PCA model) might take a long time to run. 


## References 
[1] Nakaya HI, Wrammert J, Lee EK, Racioppi L, Marie-Kunze S, et al. Systems biology of vaccination for seasonal influenza in humans. Nature immunology. 2011;12(8):786–95. https://doi.org/10.1038/ni.2067. <br/>
[2] Van Deun K, Crompvoets EAV, Ceulemans E. Obtaining insights from high-dimensional data: Sparse Principal Covariates Regression. BMC Bioinformatics. 2018; 19(104). https://doi.org/10.1186/s12859-018-2114-5  
