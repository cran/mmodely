#k number of parameters estimated including intercept
correct.AIC <- function(AIC,K,n){AIC + (2*K*(K+1))/(n-K-1)} 
# from Symonds 2011 (Eq 3)

# see Symonds and Mousalli 2011 (generalized from AIC)
weight.IC <- function(IC){
  IC[is.na(IC)|is.infinite(IC)] <- 10^10 # some impossibly high number
  best.IC <- min(IC, na.rm=TRUE)
  ICdelta <- IC - best.IC
  ICweights <- exp(-.5*ICdelta)/ do.call(sum, as.list(exp(-.5*ICdelta)))
  ICweights
}


