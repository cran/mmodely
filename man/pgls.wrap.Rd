\name{pgls.wrap}
\alias{pgls.wrap}
\title{A Wrapper for PGLS model}
\description{
Print the results of an unfit PGLS model 
}
\usage{
pgls.wrap(cd,f,b,l,k,d,all.vars=names(cd$data)[-1], 
           model.no=NA, mtx.out=NA, write=TRUE,print=FALSE)
}
\arguments{
\item{cd}{a 'comparative data' object, here created by 'comp.data(phylo, df, gs.clmn)'}
\item{f}{the model formula (as acharacter string)}
\item{b}{a list of vectors of upper and lower bounds for kappa, lambda, and delta}
\item{l}{the fixed or 'ML' value for lambda}
\item{k}{the fixed or 'ML' value for kappa}
\item{d}{the fixed or 'ML' value for delta}
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

model <- 'OC ~ mass.Kg + group.size'

pgls.wrap(cd=comp,f=model,b=list(kappa=c(.3,3),lambda=c(.3,3),delta=c(.3,3)), 
           l=1,k=1,d=1,all.vars=names(cd.obj$data)[-1])

}
