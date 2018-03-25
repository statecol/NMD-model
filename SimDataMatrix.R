rm(list=ls())  # Careful! This clears all of R's memory!

set.seed(1234)

S  =100
N = rnbinom(S, size=1, prob=0.01)
quad = 50
sim.mat = matrix(NA, nrow = S, ncol = quad)
for(i in 1:S){
  sim.mat[i,] = (rmultinom(1, N[i], prob = rep(1/quad, quad)))
}

Hypothetical.data.mat = sim.mat

save(Hypothetical.data.mat, file="/Users/shen/CloudStation_YH_TJ/1046-aggregation-scale-independence/R-GitHub/Hypothetical.data.mat.rdata")
                       