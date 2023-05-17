
#####################
## MODEL SELECTION ##
#####################

## a publication quality version of plot.pgls.iters (a plot of just R2 vs AICc)
plot.pgls.R2AIC <- function(x, bests=bestBy(x, by=c('n','q','qXn','rwGsm')[4], best=c('AICc','R2.adj')[1],inverse=c(FALSE,TRUE)[1]),  bcl=rgb(1,1,1,maxColorValue=3,alpha=1), nx=2,model.as.title='',...){
 lab.line <- par()$mgp[1]  #x is a PGLSi object
 #range.AICs <- x$AICc   < max(bests$AICc) 
 #range.R2s  <- x$R2.adj > min(bests$R2.adj)#[range.R2s&range.AICs,]
 cex.n <- x$n/max(x$n) * nx;   
 xlims <- c(min(x$R2.adj[is.finite(bests$R2.adj)]),max(x$R2.adj[is.finite(x$R2.adj)]))
 ylims <- c(min(x$AICc[is.finite(x$AICc)]),max(x$AICc[is.finite(x$AICc)]))
 with(bests,                                          plot(R2.adj, AICc, pch='', main=model.as.title, cex.main=.4, xlim=xlims, ylim=ylims, ylab='', xlab='', ...));   
 with(x,                                            points(R2.adj, AICc, pch=21, bg=bcl, col='white', bty='n', cex=cex.n))
 with(bests,                                          text(R2.adj, AICc, model.no)) ## opt. prepend "a" "p" or "ap" to these numbers
 par(xpd=FALSE)
 abline(v=0, lty=2)
 mtext(side=1, line=lab.line, text=expression(adjusted(R^2)))
 mtext(side=2, line=lab.line, text='AICc')
}


###################################
## SUBSET-SELECTED MODELS' COEFS ##
###################################

sparge.modsel <- function(PC, jit.f=1, R2x=3, nx=2, n.max=max(unlist(PC$n)), zeroline=TRUE, add=FALSE, pd=0, pvs=names(PC$coefs), pvlabs=NULL, xlim=range(unlist(PC$coefs)), MA=NULL, ap=8, ac=1, ax=nx, ...){
   
   q <- length(pvs)
   for(yi in 1:q){
    xs <- PC$coefs[rev(pvs)][[yi]]
    ys <- - jitter(rep(yi,length(xs)), factor=jit.f, amount=1/10) + pd
    lwd <- PC$R2a[rev(pvs)][[yi]] * R2x;   
    cex <- PC$n[rev(pvs)][[yi]]/n.max * nx;   
    if(yi==1 & !add){
     plot(xs, ys, ylim=-c(1,q), ylab='', cex=cex, xlim=xlim, yaxt='n', pch=21, lwd=lwd, ...) 
    }else{
     points(xs, ys, cex=cex, lwd=lwd, pch=21, ...)
    }
    if(!is.null(MA)){points(x=MA[rev(pvs)[yi]], y=-yi+pd, pch=ap, cex=ax, col=ac)} # add model averages  #FIX HARD CODING OF bg
   }
   if(!is.null(pvlabs)) pvlabs <- rev(pvlabs) else pvlabs <- rev(pvs)
   axis(2, at=-1:-q, labels=pvlabs, las=1)

 if(zeroline){abline(v=0, lty=2)}
}


################################
## OTHER USEFUL GENERAL PLOTS ##
################################

# split plot to investigate confounding between modeled variables (inspired by lattice's multi-panel lattice graph)
plot.confound.grid <- function(x,Y='y',X='x',confounder='z',breaks=3, ...){
 oldpar <- par(no.readonly=TRUE);
 par(mfrow=c(1,breaks))
 on.exit(par(oldpar))
 if(length(breaks)==1)
  breaks   <- c(quantile(  x[,confounder], probs = seq(0, 1, by = 1/breaks),na.rm=TRUE))
  ecs   <- split(x,cut(  x[,confounder], breaks=breaks));
 for(i in 1:length(ecs)){ 
  with(ecs[[i]],plot(get(X),get(Y)), ...)
  with(ecs[[i]],abline(lm(get(Y)~get(X))))
 }
}


# plot with a linear model trendline (and its p-value reported as text overlay on the plot)
plot.xy.ab.p <- function(x, x.var, y.var, fit.line=TRUE, p.value=TRUE, slope=TRUE, p.col='red', plot.labels=TRUE, verbose=TRUE, ...){
  if(verbose) print(paste('x=',x.var,'; y=',y.var))
  pch<-par()$pch
  if(plot.labels!=FALSE & !is.na(plot.labels) & !is.null(plot.labels)){
    if(length(plot.labels)>1){
      if(length(plot.labels) == nrow(x)){
       label.vector <- plot.labels
      }else{
       stop("plot label vector length does not match rows in x")
      }
    }else{
     if(is.character(plot.labels)){
       if(!plot.labels %in% names(x)) stop('cannot find plot.labels in (col)names of x')
       label.vector <- with(x,get(plot.labels))
     }else{
       label.vector <- rownames(x); message("plot.labels=TRUE so using rownames(x)")
     }
    } 
    pch <- ''      
  }
  with(x, plot(get(x.var), get(y.var), xlab=x.var, ylab=y.var, pch=pch, ...))
  if(!is.null(label.vector) & pch=='') 
  with(x, text(get(x.var), get(y.var), labels=label.vector, ...))

  model.formula <- as.formula(paste(y.var,'~',x.var))
  linear.model <- lm(model.formula,data=x)
  p <- round(summary(linear.model)$coefficients[2,4],3)
  B <- round(summary(linear.model)$coefficients[2,1],3)
  B.sign <- ifelse(test=B>0, yes="+", no="-")
  if(p.value){ print(paste(B.sign,"p-value=",p), digits=3)
               with(x,text(mean(get(x.var),na.rm=TRUE),
                 mean(get(y.var),na.rm=TRUE),labels=paste('p=',p),col=p.col)) } 
  if(slope){ print("slope=",B, digits=3)
               with(x,text(mean(get(x.var),na.rm=TRUE),
                 mean(get(y.var),na.rm=TRUE),labels=paste('slope=',B),col=p.col)) } 
  if(fit.line) with(x, abline(linear.model, col=p.col))
}
