\name{average.fit.models}
\alias{average.fit.models}
\alias{variable.importance}
\title{Calculate a weighted average of pglm}
\description{
These function takes the output of pgls.iter and uses its list of objects { model fits, optimzations (e.g.AICc)} and performs a weighted average on the ctoefficients estimated in the former by weighting by the latter. The parameters can also optionally be converted to binary by specifying "binary=FALSE" or just running the alias wraper function for assessing evidence of variable importance (Burnham & Anderson 2000).
}
\usage{
average.fit.models(vars, fits, optims,weight='AICw',
         by=c('n','q','nXq','rwGsm')[1], round.digits=5, binary=FALSE, standardize=FALSE)
variable.importance(vars, fits, optims,weight='AICw',
         by=c('n','q','nXq','rwGsm')[1], round.digits=5)
}
\arguments{
\item{vars}{variable names of model}
\item{fits}{a list of PGLS model fits}
\item{optims}{a list of PGLS optimization paramters (should include "AICw")}
\item{weight}{a column name in the optims that specifies the weights to be used in the average}
\item{by}{unique identifier used to group sub-datasets for reporting (defaults to n) }
\item{round.digits}{the tnumber of decimal places of the resultant mean to ouput}
\item{binary}{converts all parameters to binary for presense or absense to calculate 'importance'}
\item{standardize}{standardize the coefficient estimates by partial standard deviations, according to Cade (2015)}
}
\value{
A vector of AICc difference weighted [AICw] averages of PGLS coefficients.
Also returns model 'selection' errors or the square root of 'uncertainties' (Burnham & Anderson 2000) 
}
\examples{

data.path <- system.file("extdata","primate-example.data.csv", package="mmodely")
data <- read.csv(data.path, row.names=1)
pvs <- names(data[3:5])
data$gn_sp <- rownames(data)

tree.path <- system.file("extdata","primate-springer.2012.tre", package="mmodely")
phyl <- ape::read.tree(tree.path)[[5]]

comp <- comp.data(phylo=phyl, df=data)

mods <- get.model.combos(predictor.vars=pvs, outcome.var='OC', min.q=2)

PGLSi <- pgls.iter(models=mods, phylo=phyl, df=data, k=1,l=1,d=1) 

average.fit.models(vars=c('mass.Kg','group.size'), fits=PGLSi$fits, optims=PGLSi$optim)
variable.importance(vars=c('mass.Kg','group.size'), fits=PGLSi$fits, optims=PGLSi$optim)

}
