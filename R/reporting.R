
.clean.report.string <- function(string){   
  for(pat in c(c('T\\(','TR\\(','TRU\\(','TRUE\\('),'\\..*\\(')){ ### LAST ONE ("\\.") IS SPECIFIC TO MY MODELS -> GENERALIZE?
    string <- gsub(x=string, pattern=pat, replacement='(') }
  return(string)
}

fit.1ln.rprt <- function(fit, method=c('std.dev','p-value')[1], decimal.places=3, name.char.len=6, print.inline=TRUE, rtrn.line=FALSE, R2AIC=TRUE,mn=''){
  ps <- summary(fit)$coefficients[-1,4]
  ts <- summary(fit)$coefficients[-1,3]
  R2adj<- round(summary(fit)$adj.r.squared,3)
  AICc <- round(correct.AIC(AIC=AIC(fit),K=fit$k, n=fit$n), 2)
  coefs <- coefficients(fit)[-1]
  if(method=='p-values'){
   sigs <- sapply(-ceiling(log(ps,10)), function(t) paste(rep('+',min(t,3)),collapse=''))
  }else{
   sigs <- sapply(floor(abs(ts)), function(t) paste(rep('+',min(t,3)),collapse=''))
  }
  sigs[coefs<0] <- gsub(x=sigs[coefs<0],pattern='\\+','-') #seems backwards...
  names(ps) <- paste(sigs,names(coefs),sep='')
  names(ps) <- substr(names(ps),1,name.char.len)
  coefs.pos <- coefs[coefs>0]
     ps.pos <-    ps[coefs>0]
  coefs.neg <- coefs[coefs<0]
     ps.neg <-    ps[coefs<0]

  ps.p.ord <- order(ps.pos); 
  ps.n.ord <- rev(order(ps.neg))
  out.vect <- round(c(ps.pos[ps.p.ord], ps.neg[ps.n.ord]),decimal.places)
  # consider cleaning up names here so that anything after a "." is removed (as well as "TRUE")

    .paste.pvals2prnds <- function(ov){ paste(names(ov),paste("(",ov,")",sep=''),sep='')}
    report.pos <- .paste.pvals2prnds(out.vect[                 1: length(ps.pos)                ]); if(length(ps.pos)==0){report.pos <-''}
    report.neg <- .paste.pvals2prnds(out.vect[(length(ps.pos)+1):(length(ps.pos)+length(ps.neg))]); if(length(ps.neg)==0){report.neg <-''}
    out.line <- c(report.pos," | ",report.neg)
    outline <- paste(.clean.report.string(out.line),collapse=' ')
    if(R2AIC){outline <- paste(outline, 'R2adj:',R2adj,'AICc:',AICc)}
  if(print.inline){
    cat(paste(mn,outline,'\n'))#print(, quote=FALSE)
  }
  if(rtrn.line)
   invisible(outline)
  else 
   invisible(out.vect)
}


######################
### PGLS REPORTING ###
######################


pgls.print <- function(pgls, all.vars=names(pgls$data$data)[-1], model.no=NA, mtx.out=NA, write=TRUE, print=FALSE){
 R2<- summary(pgls)$r.squared
 R2adj<- summary(pgls)$adj.r.squared
 AIC <- pgls$aic   # AIC(pgls);
 AICc <- pgls$aicc # correct.AIC(AIC, K=pgls$k, n=pgls$n)
 #print(anova(pgls)); 
 if(!is.na(mtx.out)){ 
   header.row.ct <- 5
   y <- length(all.vars)+ header.row.ct    # BUG WITH INTERACTION TERMS THAT MIGHT NOT MATCH DATA WIDTH
   mtx <- matrix(nrow=y, ncol=7)

   mtx[1,1] <- 'model #';  mtx[1,3] <- model.no;     
   mtx[2:y,2] <- ''; mtx[2:y,4] <-'(';mtx[2:y,6] <- ')'
   mtx[2,1] <- 'data/model size';  mtx[2,3] <- paste('n=',pgls$n,sep='');     mtx[2,5] <- paste('k=',count.mod.vars(pgls$formula),sep='')
   mtx[3,1] <- 'R2 (adj)';          mtx[3,3] <- round(R2,3);  mtx[3,5] <- round(R2adj,3)
   mtx[4,1] <- 'AIC (AICc)';        mtx[4,3] <- round(AIC,2); mtx[4,5] <-round(AICc,2)
   mtx[5,1] <- paste(sapply(1:3,function(i) paste(names(pgls$param)[i],round(pgls$param[i],2),sep='=')),collapse=';')

   mtx[6:y,1] <- all.vars #
   wch <- unlist(sapply(names(pgls$model$coef), function(v) which(all.vars==sub(x=v,'TRUE',''))))+header.row.ct
   mtx[wch,3] <- round(summary(pgls)$coefficients[-1,'Estimate'],3)
   Ps <- summary(pgls)$coefficients[-1,'Pr(>|t|)']
   mtx[wch,5] <- round(Ps,5)
   mtx[wch,7] <- c(NA,'***')[(Ps<.001)+1];na <- which(is.na(mtx[wch,7]))
   mtx[wch[na],7] <- c(NA,'**') [(Ps[na]<=.01)+1];na <- which(is.na(mtx[wch,7]))
   mtx[wch[na],7] <- c(NA,'*')  [(Ps[na]<=.05)+1];na <- which(is.na(mtx[wch,7]))
   mtx[wch[na],7] <- c(NA,"'.")  [(Ps[na]<=.1)+1];na <- which(is.na(mtx[wch,7]))
   mtx[wch[na],7] <- c(NA,'_')  [(Ps[na]<=.5)+1];na <- which(is.na(mtx[wch,7]))

   mtx[is.na(mtx[,5]),c(4,6)] <- ''
   mtx[is.na(mtx)] <- ''

   if(print) print(mtx)
   if(write) write.csv(x=mtx,file=mtx.out)
 }
}

pgls.wrap <- function(cd,f,b,l,k,d,all.vars=names(cd$data)[-1], model.no=NA, mtx.out=NA, write=TRUE,print=FALSE){
 pgls <- pgls(formula=formula(f), lambda=l, kappa=k, delta=d, bounds=b,data=cd);
 pgls.print(pgls, all.vars=all.vars, model.no=model.no, mtx.out=mtx.out, write=write, print=print)
}

pgls.report <- function(cd, f=formula('y~x'), l=1,k=1,d=1,bounds=list(lambda=c(.2,1),kappa=c(.2,2.7),delta=c(.2,2.7)), 
                         anova=FALSE, mod.no='NA', out='pgls.output-temp',QC.plot=FALSE){ #cd: "comparative data'
 oldpar <- par(no.readonly=TRUE); on.exit(par(oldpar))
 write.to.disk <- out!=''
 pgls <- pgls(data=cd, formula=f, lambda=l, kappa=k, delta=d, bounds=bounds);
 if(QC.plot){ par(mfrow=c(2,2));plot(pgls)}
 print(summary(pgls)); print(paste("AIC =",round(AIC(pgls),1))); if(anova) print(anova(pgls)); 
 pgls.print(pgls, model.no=mod.no, mtx.out=paste(out,'-',mod.no,'.csv',sep=''), write=write.to.disk)
 sorted.p.values <- fit.1ln.rprt(pgls, R2AIC=FALSE)
 return(pgls)
}

