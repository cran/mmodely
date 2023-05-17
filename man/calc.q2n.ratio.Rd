\name{calc.q2n.ratio}
\alias{calc.q2n.ratio}
\title{Calculate the ratio of fit predictor variables to sample size}
\description{
The one in ten rule of thumb for model fitting suggest at least 10 fold as many data as parametes fit. This function allows for easily calculating that ratio on model selected PGLS fits.
}
\usage{
calc.q2n.ratio(coefs)
}
\arguments{
\item{coefs}{a list of coefficients extracted from fit PGLS models}
}
\value{
the ratio of q to n (on average for all extracted fit models) 
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

calc.q2n.ratio(coefs.objs)

}
