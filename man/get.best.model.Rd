\name{get.best.model}
\alias{get.best.model}
\title{Get the best model from list of PGLS model fits}
\description{
Get the outcome variable from the front of a model formula string. Used as part of 'get.mod.clmns' to be passed to 'comp.data'
}
\usage{
get.best.model(PGLSi, by=c('AICc','R2.adj')[1])
}
\arguments{
\item{PGLSi}{a list of PGLS iter objects, each of which is list including: fitted PGLS model, a optim table, and a tree-transformation parameter table}
\item{by}{metric to use in searching for the best model}
}
\value{
a line corresponding to the "best" model from the PGLSi "optim" table
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

a.PGLS <- get.best.model(PGLSi, by=c('R2.adj','AICc')[1])

}
