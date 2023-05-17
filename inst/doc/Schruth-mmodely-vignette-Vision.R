### R code from vignette source 'Schruth-mmodely-vignette-Vision.Rnw'
### Encoding: NA

###################################################
### code chunk number 1: foo
###################################################
opt.old <- options(keep.source = TRUE, width = 95)
foo <- packageDescription("mmodely")


###################################################
### code chunk number 2: downloadpackages
###################################################
# wget https://cran.r-project.org/src/contrib/Archive/caroline/caroline_0.8.0.tar.gz
# wget https://cran.r-project.org/src/contrib/Archive/caper/caper_0.5.tar.gz
# wget https://cran.r-project.org/src/contrib/Archive/ape/ape_3.0-5.tar.gz
# R CMD INSTALL caroline_0.8.0.tar.gz 
# R CMD INSTALL caper_0.5.tar.gz 
# R CMD INSTALL ape_3.0-5.tar.gz


###################################################
### code chunk number 3: loadlib
###################################################
library(caper)
library(mmodely)


###################################################
### code chunk number 4: read
###################################################
data.path <- system.file("extdata","primate-example.data.csv", package="mmodely")
data <- read.csv(data.path, row.names=1)
data$gn_sp <- rownames(data)

# merge data sets here if applicable

tree.path <- system.file("extdata","primate-springer.2012.tre", package="mmodely")
phyl <- ape::read.tree(tree.path)[[5]]
#5. RAxML phylogram based on the 61199 bp concatenation of 69 nuclear and ten mitochondrial genes. 

phyl <- trim.phylo(phylo=phyl, gs.vect=data$gn_sp) # prune unused nodes and branches

comp <- comp.data(phylo=phyl, df=data)


###################################################
### code chunk number 5: ti
###################################################
model <- as.formula('OC ~ mass.Kg + group.size')
fit <- caper::pgls(formula=model, data=comp)
summary(fit)

pgls.report(comp, f=model, anova=TRUE, QC.plot=TRUE)


###################################################
### code chunk number 6: iter
###################################################
pv0 <- c("mass.Kg","group.size","arboreal","nocturnal") #"swing.pct"
est.mods <- get.model.combos(predictor.vars=pv0, outcome.var='OC', min.q=2)

ps <- get.phylo.stats(phylo=phyl, data=data, trait.clmn='OC'); 
lambda <- ps$lambda$lambda ; print(lambda)
PGLSi <- pgls.iter(models=est.mods, phylo=phyl, df=data, l=lambda, k='ML', d='ML')

pgls.iter.stats(PGLSi) # check run, especially to see how few sub-datasets exist


###################################################
### code chunk number 7: ttavg
###################################################
tt.avgs <- apply(PGLSi$params, 2, mean, na.rm=TRUE) # tree transformation averages
print(tt.avgs)


###################################################
### code chunk number 8: fixiter
###################################################

pvs <- c("mass.Kg","group.size","infant.carry","arboreal","DPL.km","nocturnal")
all.mods <- get.model.combos(predictor.vars=pvs, outcome.var='OC', min.q=2)

data <- subset(data,!grepl(rownames(data),pattern='gorilla')) # remove an OC measurement outlier

# randomly sprinkle in some missing values (for more interesting for model selection)
missing.value.ct <- 1
for(pv in pv0){ data[sample(x=1:nrow(data),size=missing.value.ct),pv] <- NA} 

PGLSi <- pgls.iter(models=all.mods, phylo=phyl, df=data, l=lambda, k=tt.avgs['k'], d=tt.avgs['d'])


###################################################
### code chunk number 9: fixiterstats
###################################################
pgls.iter.stats(PGLSi) 


###################################################
### code chunk number 10: modavg
###################################################
w.means.pds <- average.fit.models(vars=pvs, fits=PGLSi$fits, optims=PGLSi$optim, by='rwGsm', standardize=TRUE)
#
apply(w.means.pds, 2, mean, na.rm=T) #average of weighted means over all sub-datasets
w.means.pds                                    # weighted means    per   sub-dataset


###################################################
### code chunk number 11: varimport
###################################################
w.import.pds <- variable.importance(vars=pvs, fits=PGLSi$fits, optims=PGLSi$optim, by='rwGsm')
#
apply(w.import.pds, 2, mean, na.rm=T) #average of weighted means over all sub-datasets
w.import.pds                                    # weighted means    per   sub-dataset


###################################################
### code chunk number 12: modsel
###################################################
select.best.models(PGLSi, using='AICc') 


###################################################
### code chunk number 13: plotiter
###################################################
plot.pgls.iters(PGLSi)


###################################################
### code chunk number 14: plotiter
###################################################
plot.pgls.iters(PGLSi)


###################################################
### code chunk number 15: getcoef
###################################################
sdevs.objs <- get.pgls.coefs(PGLSi$fits, est='t value')
coefs.objs <- get.pgls.coefs(PGLSi$fits, est='Estimate')


###################################################
### code chunk number 16: lnrpt
###################################################
report.vect <- sapply(1:length(PGLSi$fits), function(i) fit.1ln.rprt(PGLSi$fits[[i]], rtrn.line=FALSE, mn=i))


###################################################
### code chunk number 17: R2AIC
###################################################
par(mar=c(5,5,3,3))
plot.pgls.R2AIC(PGLSi$optim)


###################################################
### code chunk number 18: R2AIC
###################################################
par(mar=c(5,5,3,3))
plot.pgls.R2AIC(PGLSi$optim)


###################################################
### code chunk number 19: coefdists
###################################################
par.old <- par(mar=c(5,8,1,4),mfrow=c(2,1))
sparge.modsel(sdevs.objs, R2x=7, xlab='t value')
sparge.modsel(coefs.objs, R2x=7, xlab='Estimate')


###################################################
### code chunk number 20: coefdists
###################################################
par.old <- par(mar=c(5,8,1,4),mfrow=c(2,1))
sparge.modsel(sdevs.objs, R2x=7, xlab='t value')
sparge.modsel(coefs.objs, R2x=7, xlab='Estimate')


###################################################
### code chunk number 21: bar
###################################################
options(opt.old)
par(par.old)


