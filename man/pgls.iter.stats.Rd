\name{pgls.iter.stats}
\alias{pgls.iter.stats}
\title{Statistics from PGLS runs}
\description{
Print (and plot) statistics from a list of PGLSs fitted models and tables of associated parameters.
}
\usage{
pgls.iter.stats(PGLSi, verbose=TRUE, plots=FALSE)
}
\arguments{
\item{PGLSi}{a list of PGLS iter objects, each of which is list including: fitted PGLS model, a optim table, and a tree-transformation parameter table}
\item{verbose}{the model formula (as acharacter string)}
\item{plots}{the fixed or 'ML' value for kappa}
}
\value{
A summary statistics on each of the objects in the PGLS list of lists
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

pgls.iter.stats(PGLSi, verbose=TRUE, plots=FALSE)

}
