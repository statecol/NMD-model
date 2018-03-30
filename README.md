# NMD-model
This R code is to fit the negative multinomial model (NMD) or negative binomial model (NBD) for modeling community-level multi-species correlated distribution patterns in ecology. 

The similarity between NMD and NBD models is that both will be identical when local species abundance dataset is merged together.
The difference between NMD and NBD models is that NMD model is suited to model nonrandom or correlated distribution of species across different local samples.    

Before using the code here, it is necessary to introduct a bit about the documents stored here:    

1, NMD-or-NBD-estimation.R file contains all the functions for implementing and fitting the NMD or NBD model. 

2, Hypothetical.data.mat.rdata is a hypothetical species-site abundance dataset in R data file format.  The data file contains 100 species across 50 sites. This dataset demonstrates an example of nonrandom distribution of species or local species abundance distribution (local SAD). 

3, bci.Xi.rdata is another R data file containing the empirical abundance of 300 tree species recorded in BCI forest plot (2005 census). This dataset demonstrates an example of regional species abundance distribution (regional SAD). 

To use the R code, it is recommended to download all the files (particularly the above three files: NMD-or-NBD-estimation.R, Hypothetical.data.mat.rdata and bci.Xi.rdata), save them into the same local directory. Open R and type the following R commands: 

 
 
 
=================================================================================

setwd("your path to the local working directory of the downloaded rdata and R files") 


source("NMD-or-NBD-estimation.R") 

load('bci.Xi.rdata') 

Xi = bci.Xi 

xmax = max(Xi) 

x = factor(Xi, levels = 1:xmax) 

f = table(x) 

print(NMD.NBD.Estimation(ini.par = c(0.2, 0.2), f=f, m=1))
 
 
 
 
load('Hypothetical.data.mat.rdata') 

print(NMD.DataMatrix.Estimation(ini.par = c(0.2, 0.2), Xi.mat=Hypothetical.data.mat, m=1)) 

