# simple wrapper function for "comparative.data" function that keeps vcv info consistent (used in pgls.print and pgls.iter)
comp.data <- function(phylo,df,gn_sp='gn_sp'){ comparative.data(phy=phylo, data=df, names.col=gn_sp, vcv=TRUE, vcv.dim=3)}#, warn.dropped=TRUE)} #noNA(df,vars)v, 

# requires caper's comp.data and pgls
pgls.iter <- function(models, phylo, df, gs.clmn='gn_sp', b=list(lambda=c(.2,1),kappa=c(.2,2.7),delta=c(.2,2.7)), l='ML',k='ML',d='ML'){ #, vcv.d=3){

 ### iterate through these possible combinations of these params for ^R2 & vAIC
 fits <- list()
 params <- matrix(rep(c(NA,NA,NA), length(models)), ncol=3);
 R2s <- R2s.adj <- AICs <- AICcs <-  rwGsm <- nv(rep(NA, length(models)) ,models) #model.lengths <-

 rownames(params)<-models; #colnames(params) <- c('l','k','d');
 for(model in models){
   cat(paste(which(models %in% model), model,'\n'))#print(, quote=FALSE)
   cd <-  comp.data(      phylo=phylo,  df=df[,c(get.mod.clmns(model),gs.clmn)], gs.clmn)#gn_sp   # use drop.na.data here?
   #cd <- comparative.data(phy  =phylo,data=df[,get.mod.clmns(model)], names.col=gs.clmn, vcv=vcv.dim>0, vcv.dim=vcv.d) 
   fits[[model]] <- pgls(formula=formula(model), lambda=l, kappa=k, delta=d, data=cd, bounds=b)  
   params[model,] <- fits[[model]]$param[names(b)]
   R2s[model] <- summary(fits[[model]])$r.squared
   R2s.adj[model] <- summary(fits[[model]])$adj.r.squared
   AICs[model] <- fits[[model]]$aic #AIC(fits[[model]])
   AICcs[model] <- fits[[model]]$aicc #correct.AIC(AIC=AIC(fits[[model]]), K=fits[[model]]$k, n=fits[[model]]$n) 
   #model.lengths[model] <- sapply(strsplit(model,split='\\+'),length)
 }

 optim.1 <- data.frame(n=sapply(fits,function(i) i$n), n=sapply(fits,function(i) i$n),     ## DOES $n include rows with MISSING DATA?
                       q=sapply(models, function(m) count.mod.vars(formula(m)))) #count the predictor variables
 optim.1$qXn <- apply(optim.1[,c('q','n')], 1, paste, collapse='X')# for all the unique combinations of n and q
 optim.1$rwGsm <-sapply(models, function(m) sum(as.numeric(sapply(strsplit(paste(rownames( ##  drop.na.data removes rows w/ MISSING DATA
                      drop.na.data(fits[[m]]$data$data[,fits[[m]]$varNames])),collapse='.'),'')[[1]], charToRaw))))# unique # for each G_sXn combo

 optim.2 <- data.frame(model.no=1:length(models), R2=R2s,R2.adj=R2s.adj,AIC=AICs,AICc=AICcs, stringsAsFactors=FALSE); 
 optim.2$AICw <- weight.AIC(AIC=optim.2$AICc);

 optim <- cbind.data.frame(optim.1,optim.2)

 row.names(optim) <- models
 colnames(params) <- substr(x=names(b),1,1) #c('l','k','d');

 return(list(fits=fits,params=params,optim=optim))

}

pgls.iter.stats <- function(PGLSi, verbose=TRUE, plots=FALSE){
 oldpar <- par(no.readonly=TRUE); 
                            # PGLSi object switched to "x" in plotting functions
 optim <- PGLSi$optim
 dataset.dims <-  apply(optim[,c('q','n','qXn','rwGsm')], 2, function(c){length(unique(c))} )
 dataset.avgs <-  apply(optim[,c('q','n')], 2, summary)
 optim.quants <- optim[,sapply(optim[,-ncol(optim)], function(x) {!is.character(x)})]
 param.stats  <-   apply(PGLSi$params, 2, mean, na.rm=TRUE)

 if(verbose){ 
  cat(paste('models:',length(PGLSi$fits)),sep='\n')
  cat('dimensions of sub-datasets:',sep='\n')
   print(dataset.dims)
   print(dataset.avgs)
  cat('tree transformation parameter averages:',sep='\n')
   print(param.stats)
  cat('distributions of optimization parameters:',sep='\n')
  print(summary(optim.quants) )
 }
 if(plots){
  par(mfrow=c(3,3));on.exit(par(oldpar))
  sapply(names(optim.quants), function(oq) hist(optim.quants[,oq], main=oq))
 }
}
#######################
### MODEL AVERAGING ###
#######################

average.fit.models <- function(vars, fits, optims, weight='AICw', by=c('n','q','nXq','rwGsm')[4], round.digits=5){
  q <- length(vars)
  f <- length(fits)
  coef.mtrx <- matrix(data=rep(NA,q*f), ncol=q,nrow=f);
  rownames(coef.mtrx)<- names(fits)
  colnames(coef.mtrx)<-vars
  for(fit.nm in names(fits)){ 
   coef.this <- coef(fits[[fit.nm]])
   names(coef.this) <- sapply(names(coef.this),function(cf) sub(x=cf,'TRUE',''))
   for(vp in vars){ 
    coef.mtrx[fit.nm,vp] <- coef.this[vp]
   }
  }
  # weighted means of coefficients based on AIC
  coef.opt <- cbind.data.frame(list(coef.mtrx,optims))
  coef.opts <- split(x=coef.opt, f=as.factor(coef.opt[,by]))
           # old version -> #round(sapply(vars, function(vi) weighted.mean(x=coef.mtrx[,vi],w=optims[,weight],na.rm=TRUE)),round.dec.pl) 
  ret <- sapply(coef.opts, function(co) round(sapply(vars, function(vi) weighted.mean(x=co[,vi],w=co[,weight],na.rm=TRUE)),round.digits))
  return(t(ret))
}

######################
## MODEL SELECTION ###
######################


select.best.models <- function(PGLSi, using=c('AICc','R2.adj','AIC','R2')[1], by=c('n','q','nXq','rwGsm')[4]){ 
     if(grepl(x=using,'AIC')){   
        ret <- bestBy(df=PGLSi$optim, by=by, best='AICc')
     }else{
        ret <- bestBy(df=PGLSi$optim, by=by, best='R2.adj', inverse=T)
     }
 return(ret)
}


plot.pgls.iters <- function(x, bests=bestBy(x$optim, by=c('n','q','qXn','rwGsm')[4], best=c('AICc','R2.adj')[1],inverse=FALSE), ...){
 oldpar <- par(no.readonly=TRUE);
 par(mfrow=c(2,2),mar=c(3,3,1.5,1.5),mgp=c(1.7,0.4,0))       #x is a PGLSi object            
 on.exit(par(oldpar))                 
 with(x$optim, plot(R2, AIC))
 with(bests, plot(R2,AIC,pch=''));with(bests, text(R2,AIC,model.no))
 with(x$optim, plot(R2.adj, AICc))
 with(bests, plot(R2.adj,AICc,pch=''));with(bests, text(R2.adj,AICc,model.no))
}



###################################
## SUBSET-SELECTED MODELS' COEFS ##
###################################

get.pgls.coefs <- function(pgls.fits, est=c("t value","Estimate","Pr(>|t|)")[1]){  
 coef.lst <- R2a.lst <- list() # nv(rep(NA,length(coef.names)),coef.names))
 var.names=unique(unlist(lapply(1:length(pgls.fits), function(i) pgls.fits[[i]]$varNames)))[-1]
  for(pgls.fit in pgls.fits){
   coef.names <- pgls.fit$varNames[-1]
   for(vi in 1:length(coef.names)){ 
    coef.lst[[coef.names[vi]]] <- c(coef.lst[[coef.names[vi]]], summary(pgls.fit)$coefficients[vi+1,est])
     R2a.lst[[coef.names[vi]]] <- c( R2a.lst[[coef.names[vi]]], summary(pgls.fit)$adj.r.squared)
    }
   }
  return(list(var.names=var.names, coefs=coef.lst, R2adj=R2a.lst))
}

