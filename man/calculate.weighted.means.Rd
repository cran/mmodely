\name{calculate.weighted.means}
\alias{calculate.weighted.means}
\title{Calculate a weighted average of pglm}
\description{
This function takes the output of pgls.iter and uses its list of objects { model fits, optimzations (e.g.AICc)} and performs a weighted average on the ctoefficients estimated in the former by weighting by the latter.
}
\usage{
calculate.weighted.means(vars, fits, optims, weight='AICw', decimal.places=5)
}
\arguments{
\item{vars}{variable names of model}
\item{fits}{a list of PGLS model fits}
\item{optims}{a list of PGLS optimization paramters (should include "AICw")}
\item{weight}{a column name in the optims that specifies the weights to be used in the average}
\item{decimal.places}{the tnumber of decimal places of the resultant mean to ouput}
}
\value{
A vector of AICc difference weighted [AICw] averages of PGLS coefficients
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

calculate.weighted.means(vars=c('mass.Kg','group.size'), fits=PGLSi$fits, optims=PGLSi$optim)

}
