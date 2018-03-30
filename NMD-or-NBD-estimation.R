library(nlme)


# The following R function cond.logf() is the log likelihood function 
#   whose antilogarithm is the likelihood function in Eq. 8b of the main text.

# There are three arguments in this R function:
#   The first one is a numeric vector of the two parameters a and b.
#       a: the aggregation parameter denoted by k in the main text
#       b: it is a reparametrized parameter by b = m*u/a 
#   The second one is a numeric vector of the species frequency counts: f1, f2, ...
#       in which f1 is the number of singletons, f2 is the number of doubletons, and so forth.
#   The third one is m: the area size of the local sample; the default was set at m=1. 
#   Note that the value of m does not influence the estimation of k and u. 

############### Log likelihood function #################
cond.logf <- function(x, f, m=1) {
  a <- x[1]
  b <- x[2]
  
  ## 1-pp: the success probability of NBD model
  pp = b/(b + m)
  
  ## Determine the indices of nonzero species frequency counts
  zz = which(f > 0)  
  
  rhoN = lgamma(zz + a) - lgamma(zz + 1) - lgamma(a) + zz * 
    log1p(-pp) + a * log1p(pp - 1) - log(1 - pp^a)
  
  ## This term is solely related to observed data, 
  ## and thus it has no effect on the MLEs of parameters.
  ## However, it needs to include this term when calculating AIC. 
  const = lgamma(sum(f)+1)-sum(lgamma(f+1))
  
  res <- -sum(f[zz] * rhoN) - const
  
  res
}


# This function is to convert a vector of species abundance data (X) to 
# a vector of species frequency counts data (f).
X.to.f = function(X){
  f = factor(X, levels = 0:max(X))
  f = table(f, exclude = 0) ## frequency counts
  dimnames(f) = NULL
  return(f)
} 


##### main function 
## This R function NMD.NBD.Estimation() is to compute the MLEs, estimated SE, and 
## 95% confidence intervals of k and u of NMD (or NBD) when applied to 
## species abundance data or species frequency counts data.
## ini.par: the initial value of parameters a and b 
##          when numerically calculating the MLEs of a and b
## f: species frequency counts data
## xi: species abundance data
## Note that users just need to input either species frequency counts (f) or 
##  species counts (xi) in this function.
NMD.NBD.Estimation = function(ini.par = c(0.2, 0.003), f=NULL, xi=NULL, m=1){
  if(is.null(f)*is.null(xi)) {
    print("Please input either species frequency counts data or species abundance data!!")
    return(NULL);
  }
  if(is.null(f)) {
    f = X.to.f(xi)
  } 
  
  ## If the aggregation parameter k attains the lower bound (1e-05 here), 
  ## NBD model is equivalent to the logseries model, while it attains 
  ## the upper limit (1000 here), NBD model can be treated as a Poisson model.
  cond.logfSol = nlminb(start = ini.par, cond.logf, f = f, m = m, 
                        lower = c(1e-05, 1e-05), upper = c(1000, 1000))
  
  ## The resulting MLEs of a and b
  cond.ahat = cond.logfSol$par[1]
  cond.bhat = cond.logfSol$par[2]
  ## the MLE of u is back transformed from the reparametrizing formula b = m*u/a 
  cond.uhat = cond.ahat*cond.bhat/m
  
  ## Calculate the observed information matrix in Eq. S1 of the supporting information 
  Fisher.Info = (fdHess(c(cond.ahat, cond.bhat), f = f, m = m, fun = cond.logf)$Hessian)
  
  ## the estimated variance of a  
  var.a = solve(Fisher.Info)[1,1]
  ## the estimated variance of b  
  var.b = solve(Fisher.Info)[2,2]
  
  ## the estimated variance of u (mean abundance) which is derived from the transformation u = ab/m
  var.u = matrix(c(cond.bhat,cond.ahat)/m, nrow=1)%*%
    solve(Fisher.Info)%*%matrix(c(cond.bhat,cond.ahat)/m, ncol=1)
  
  ## Creating an R dataframe to store the computing results
  output = NULL
  output=rbind(output, c(round(cond.ahat,5), round(sqrt(var.a),5), 
                         round(max(cond.ahat-1.96*sqrt(var.a),0),5), round(cond.ahat+1.96*sqrt(var.a),5)))
  output=rbind(output, c(round(cond.uhat,5), round(sqrt(var.u),5), 
                         round(max(cond.uhat-1.96*sqrt(var.u),0),5), round(cond.uhat+1.96*sqrt(var.u),5)))
  
  output = data.frame(output)
  output = cbind(c("k","u"), output)

  colnames(output)=
    c("Parameter","Estimate","Estimated SE","95% lower limit","95% upper limit")
  
  output
} ### end of R function NMD.NBD.Estimation()     





# The following R function log.MD() is the log likelihood function 
#   whose antilogarithm is the likelihood function on the left part of the second equality of 
#   Eq. 6 in the main text.
# The argument x stands for abundances over sampled quadrats for a single species

#### Log Multinomial likelihood ### 
log.MD = function(x){
  N =sum(x)
  ## the number of quadrats in the sample
  q=length(x)
  
  val = lgamma(N+1)-sum(lgamma(x+1))-N*log(q)
  val
}



# The following R function NMD.logf() is the log likelihood function 
#   whose antilogarithm is the likelihood function in 
#   Eq. 6 in the main text.
# There are three arguments in the function:
#   The first one (x) is a numeric vector of the two parameters a and b.
#       a: the aggregation parameter denoted by k in the main text
#       b: it is a reparametrized parameter by b = m*u/a 
#   The second one (Xi.mat) is a numeric matrix of abundances, in which, 
#       the row index represents species label and the column index represents 
#       the quadrat index.
#   The third one (m) is the quadrat size; the default was set at m=1.
#   Note that the value of m does not influence the estimation of k and u. 
############### NMD likelihood #################
NMD.logf <- function(x, Xi.mat, m=1) {
  a <- x[1]
  b <- x[2]
  
  ## 1-pp: the success probability of NBD model
  pp = b/(b + m)

  ## the left part of Eq. 6 of the main text
  log.mult.prob = sum(apply(Xi.mat,1,log.MD))
  
  # Calculate species abundance data from the data matrix "Xi.mat"
  Xi = apply(Xi.mat,1,sum)
  # the species frequency counts: f1, f2, ...
  f = X.to.f(Xi)
  
  ## the right part of Eq. 6 of the main text when conditioned on
  ## the observed number of species
  zz = which(f > 0)  
  rhoN = lgamma(zz + a) - lgamma(zz + 1) - lgamma(a) + zz * 
    log1p(-pp) + a * log1p(pp - 1) - log(1 - pp^a)

  res <- -sum(f[zz] * rhoN)-log.mult.prob
  
  res
}





##### main function 
##  This R function NMD.DataMatrix.Estimation() is to compute the MLEs, estimated SE, and 
##  95% confidence intervals of k and u of NMD when applied to a data matrix of abundances, 
##  in which, the row index represents species label and the column index represents 
##       quadrat index.
##  There are three arguments in the function:
##  ini.par: the initial value of parameters a and b 
##          when numerically calculating the MLEs of a and b
##  The second one (Xi.mat) is a numeric matrix of abundances as introduced above.
#   The third one (m) is the quadrat size; the default was set at m=1.
NMD.DataMatrix.Estimation = function(ini.par = c(0.2, 0.003), Xi.mat, m=1){

  NMD.logfSol = nlminb(ini.par, NMD.logf, Xi.mat=Xi.mat, m = m, 
                       lower = c(1e-05, 1e-05), upper = c(1000, 1000))
    
  ## The resulting MLEs of a and b
  cond.ahat = NMD.logfSol$par[1]
  cond.bhat = NMD.logfSol$par[2]
  cond.uhat = cond.ahat*cond.bhat/m
  
  ## Calculate the observed information matrix 
  Fisher.Info = (fdHess(c(cond.ahat, cond.bhat), Xi.mat=Xi.mat, m = m, fun = NMD.logf)$Hessian)

  ## the estimated variance of a  
  var.a = solve(Fisher.Info)[1,1]
  ## the estimated variance of b  
  var.b = solve(Fisher.Info)[2,2]
  
  ## the estimated variance of u (mean abundance) which is derived from the transformation u = ab/m
  var.u = matrix(c(cond.bhat,cond.ahat)/m, nrow=1)%*%
    solve(Fisher.Info)%*%matrix(c(cond.bhat,cond.ahat)/m, ncol=1)
  
  ## Creating an R dataframe to store the computing results
  output = NULL
  output=rbind(output, c(round(cond.ahat,5), round(sqrt(var.a),5), 
                         round(max(cond.ahat-1.96*sqrt(var.a),0),5), round(cond.ahat+1.96*sqrt(var.a),5)))
  output=rbind(output, c(round(cond.uhat,5), round(sqrt(var.u),5), 
                         round(max(cond.uhat-1.96*sqrt(var.u),0),5), round(cond.uhat+1.96*sqrt(var.u),5)))
  
  output = data.frame(output)
  output = cbind(c("k","u"), output)
  
  colnames(output)=
    c("Parameter","Estimate","Estimated SE","95% lower limit","95% upper limit")
  
  output
} ### end of R function NMD.DataMatrix.Estimation()     
