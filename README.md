# NMD-model
This R code is to fit the negative multinomial model (NMD) or negative binomial model (NBD) for modeling community-level multi-species correlated distribution patterns in ecology. 

The similarity between NMD and NBD models is that both will be identical when local species abundance dataset is merged together.
The difference between NMD and NBD models is that NMD model is suited to model nonrandom or correlated distribution of species across different local samples.    

Before using the code here, it is necessary to introduct a bit about the documents stored here:    

1, NMD-or-NBD-estimation.R file contains all the functions for implementing and fitting the NMD or NBD model.
2, Hypothetical.data.mat.rdata is a hypothetical species-site abundance dataset.  The data file contains 100 species across 50 local sites.
