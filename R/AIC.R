#k number of parameters estimated including intercept
correct.AIC <- function(AIC,K,n){AIC + (2*K*(K+1))/(n-K-1)} 
# from Symonds 2011 (Eq 3)

# see Symonds and Mousalli 2011
weight.AIC <- function(AIC){
  AIC[is.na(AIC)|is.infinite(AIC)] <- 10^10 # some impossibly high number
  best.AIC <- min(AIC, na.rm=TRUE)
  AICdelta <- AIC - best.AIC
  AICweights <- exp(-.5*AICdelta)/ do.call(sum, as.list(exp(-.5*AICdelta)))
  AICweights
}


