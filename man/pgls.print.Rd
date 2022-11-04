\name{pgls.print}
\alias{pgls.print}
\title{Print the results of a PGLS model fit}
\description{
Print the results of a PGLS model fit
}
\usage{
pgls.print(pgls, all.vars=names(pgls$data$data)[-1], 
           model.no=NA, mtx.out=NA, write=TRUE, print=FALSE)
}
\arguments{
\item{pgls}{a fit PGLS model}
\item{all.vars}{the names of all the variables to be reported}
\item{model.no}{the model number (can be the order that models were run}
\item{mtx.out}{should a matrix of the tabular summary results be returned}
\item{write}{should the matrix of summary results be written to disk?}
\item{print}{should the matrix of summary results be printed to screen?}
}
\value{
A matrix of summary results of a fit PGLS model
}
\examples{

data.path <- system.file("extdata","primate-example.data.csv", package="mmodely")
data <- read.csv(data.path, row.names=1)
data$gn_sp <- rownames(data)

tree.path <- system.file("extdata","primate-springer.2012.tre", package="mmodely")

#5. RAxML phylogram based on the 61199 bp concatenation of 69 nuclear and ten mitochondrial genes.
phyl <- ape::read.tree(tree.path)[[5]]

phyl <- trim.phylo(phylo=phyl, gs.vect=data$gn_sp)

comp <- comp.data(phylo=phyl, df=data)

a.PGLS <- caper::pgls(formula('OC~mass.Kg + DPL.km'),  data=comp)

pgls.print(a.PGLS, all.vars=names(a.PGLS$data$data)[-1],
           model.no=NA, mtx.out='', write=FALSE, print=FALSE)


}
