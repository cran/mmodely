# Variance Inflation Factor (Cade 2015) rewriten for pgls
.VIF <- function(FM, FMs=NULL){
  d <- FM$data$data
  ov <- names(d)[1]
  pvs <- names(d)[-1]
  if(length(pvs==1)){
   return(1)
  }else{
  d <- cbind.data.frame(d, gn_sp=rownames(d))
  R2s <- nv(rep(NA,length(pvs)), pvs)
  for(pv in pvs){  # admitedly much slower than might be using an existing vcov ...
   f <- paste(paste(pvs[pvs!=pv], collapse='+'),'~',pv,sep='')
   if(is.null(FMs)){
    PGLSi <- pgls(formula=as.formula(f), data=comp.data(phylo=FM$data$phy, df=d), lambda=FM$param['lambda'], kappa=FM$param['kappa'], delta=FM$param['delta'], bounds=FM$bounds) 
   }else{
    if(!(f %in% names(FMs))){stop(paste(f,"does not appear in named list of fit models"))}
    PGLSi <- FMs[[as.character(f)]]
   }
   R2s[pv] <- summary(PGLSi)$r.squared
  }
  return(1/(1-R2s))
  }
}

# verbatim reproduction of MuMIn::.partialsd()); originally written by Kamil Barton 
.part.sd <- function(sd, vif, n, p) {
	sd[-1] * sqrt(1 / vif) * sqrt((n - 1) / (n - p))  # equation devised by Cade 2015    # sum of Akaike weights for each param instead?
}

# an implementation of standardization proposed by Cade 2015
.partial.sd <- function(FM, FMs=NULL){
 # a simplied version of MuMIn::std.coef() originally written by Kamil Barton 
 sd.p <- .part.sd(apply(FM$data$data, 2, sd),  vif=.VIF(FM, FMs), n=FM$n, p=length(coef(FM))-1) # see 'MuMIn' by Kamil Barton
}

# function that actually performs the multiplication of the ratio to the pre-existing coefficients
.std.coefs.by.partial.sd <-function(f, FMs){   
  sd.part <- .partial.sd(FMs[[f]]                                   , FMs); 
  sd.full <- .partial.sd(FMs[[which.max(sapply(names(FMs), nchar))]], FMs) #the first model should be the full model
   B.part <-        coef(FMs[[f]])[-1] * (sd.part/sd.full[names(sd.part)])   #remove the intercept from FM and divide part by (those parts in) full 


  return(B.part)
}


