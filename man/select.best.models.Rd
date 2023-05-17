\name{select.best.models}
\alias{select.best.models}
\title{Get the best model from list of PGLS model fits}
\description{
Get the outcome variable from the front of a model formula string. Used as part of 'get.mod.clmns' to be passed to 'comp.data'
}
\usage{
select.best.models(PGLSi, using=c('AICc','R2.adj','AIC','R2')[1],
                             by=c('n','q','nXq','rwGsm')[1])
}
\arguments{
\item{PGLSi}{a list of PGLS iter objects, each of which is list including: fitted PGLS model, a optim table, and a tree-transformation parameter table}
\item{using}{performance metric to use in searching for the best model}
\item{by}{unique identifier used to group sub-datasets for reporting (defaults to n) }
}
\value{
a line corresponding to the "best" models from the PGLSi "optim" table
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

a.PGLS <- select.best.models(PGLSi, by=c('R2.adj','AICc')[1])

}
