\name{pgls.iter}
\alias{pgls.iter}
\title{Iterate through PGLS estimations}
\description{
This function takes phylogenetic tree and a list of (all possible) combinations of variables as a vector of model strings and estimates PGLS fits based on the bounds or tree parameters provided seperately.  
}
\usage{
pgls.iter(models, phylo, df,  gs.clmn='gn_sp', 
          b=list(lambda=c(.2,1),kappa=c(.2,2.7),delta=c(.2,2.7)),l='ML', k='ML',d='ML') 
}
\arguments{
\item{models}{a vector of all possible model formulas (as character strings)}
\item{phylo}{a phylogenetic tree}
\item{df}{the name of the column used to specify 'Genus_species'}
\item{gs.clmn}{the name of the column containing the 'Genus_species' vector}
\item{b}{a list of vectors of upper and lower bounds for kappa, lambda, and delta}
\item{k}{the fixed or 'ML' value for kappa}
\item{l}{the fixed or 'ML' value for lambda}
\item{d}{the fixed or 'ML' value for delta}
}
\value{
a list of fit PGLS regression models plus 'optim' and 'param' support tables 
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

}
