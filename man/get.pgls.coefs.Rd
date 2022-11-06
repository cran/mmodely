\name{get.pgls.coefs}
\alias{get.pgls.coefs}
\title{Get coeficients from a list of PGLS model-fits (from each selected subset)}
\description{
Post PGLS model selection, the list of all possible PGLS model fits can be subset and passed to this function, which harvests out the coefficients or t-values for each model into bins for the coefficients
}
\usage{
get.pgls.coefs(pgls.fits, est=c("t value","Estimate","Pr(>|t|)")[1])
}
\arguments{
\item{pgls.fits}{a list of PGLS models output from 'pgls' or 'pgls.report' }
\item{est}{a character string indicating if Estimate or t value should be used as data points in the plot, default is 'Estimate'}
}
\value{
A list of PGLS coeficients (lists of estimates and t-values) organized by coeficient-named bins
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

coefs.objs <- get.pgls.coefs(PGLSi$fits, est='Estimate')

}
