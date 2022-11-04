### R code from vignette source 'Schruth-mmodely-vignette-Vocal.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: foo
###################################################
opt.old <- options(keep.source = TRUE, width = 60)
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

#multiply two vocalization metrics together to create "vocal complexity"
data$VC <- apply(data[,c('syllables_max','rhythm_max')], 1, prod)  

# merge data sets here if applicable

tree.path <- system.file("extdata","primate-springer.2012.tre", package="mmodely")
phyl <- ape::read.tree(tree.path)[[5]]
#5. RAxML phylogram based on the 61199 bp concatenation of 69 nuclear and ten mitochondrial genes. 

phyl <- trim.phylo(phylo=phyl, gs.vect=data$gn_sp) # prune unused nodes and branches

comp <- comp.data(phylo=phyl, df=data)


###################################################
### code chunk number 5: ti
###################################################
model <- as.formula('VC ~ mass.Kg + group.size')
fit <- caper::pgls(formula=model, data=comp)
summary(fit)

pgls.report(comp, f=model, anova=TRUE, QC.plot=TRUE)


###################################################
### code chunk number 6: iter
###################################################
pv0 <- c("mass.Kg","arboreal","home.range","monogamy") #"swing.pct"
est.mods <- get.model.combos(predictor.vars=pv0, outcome.var='VC', min.q=2)

ps <- get.phylo.stats(phylo=phyl, data=data, trait.clmn='VC'); 
lambda <- ps$lambda$lambda ; print(lambda)
PGLSi <- pgls.iter(models=est.mods, phylo=phyl, df=data, l=lambda, k='ML', d='ML')


###################################################
### code chunk number 7: modavg
###################################################
tt.avgs <- apply(PGLSi$params, 2, mean, na.rm=TRUE) # tree transformation averages
print(tt.avgs)

pvs <- c("mass.Kg","group.size","arboreal","monogamy","leap.pct","swing.pct")
all.mods <- get.model.combos(predictor.vars=pvs, outcome.var='VC', min.q=2)


# randomly sprinkle in some missing values (to keep things interesting for model selection)
missing.value.ct <- 1
for(pv in pvs){ data[sample(x=1:nrow(data),size=missing.value.ct),pv] <- NA} 

PGLSi <- pgls.iter(models=all.mods, phylo=phyl, df=data, l=lambda, k=tt.avgs['k'], d=tt.avgs['d'])


###################################################
### code chunk number 8: modavg
###################################################
calculate.weighted.means(vars=pvs, fits=PGLSi$fits, optims=PGLSi$optim)


###################################################
### code chunk number 9: plotiter
###################################################
plot.pgls.iters(PGLSi)


###################################################
### code chunk number 10: plotiter
###################################################
plot.pgls.iters(PGLSi)


###################################################
### code chunk number 11: getcoef
###################################################
sdevs.objs <- get.pgls.coefs(PGLSi$fits, est='t value')
coefs.objs <- get.pgls.coefs(PGLSi$fits, est='Estimate')


###################################################
### code chunk number 12: lnrpt
###################################################
report.vect <- sapply(1:length(PGLSi$fits), function(i) fit.1ln.rprt(PGLSi$fits[[i]], rtrn.line=FALSE, mn=i))


###################################################
### code chunk number 13: R2AIC
###################################################
plot.pgls.R2AIC(PGLSi$optim)


###################################################
### code chunk number 14: R2AIC
###################################################
plot.pgls.R2AIC(PGLSi$optim)


###################################################
### code chunk number 15: coefdists
###################################################
par.old <- par(mar=c(5,8,1,1),mfrow=c(2,1))
modsel.distro.dots(sdevs.objs, R2x=7, xlab='t value')
modsel.distro.dots(coefs.objs, R2x=7, xlab='Estimate')


###################################################
### code chunk number 16: coefdists
###################################################
par.old <- par(mar=c(5,8,1,1),mfrow=c(2,1))
modsel.distro.dots(sdevs.objs, R2x=7, xlab='t value')
modsel.distro.dots(coefs.objs, R2x=7, xlab='Estimate')


###################################################
### code chunk number 17: bar
###################################################
options(opt.old)
par(par.old)


