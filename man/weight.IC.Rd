\name{weight.IC}
\alias{weight.IC}
\title{Get IC weights }
\description{
An implementation of IC weighting that first calulates the difference in IC values by subtracting all values from the lowest IC value. Second, the changes are expoentiated divided by a sum of the same and exponentiated yet again.
}
\usage{
weight.IC(IC)
}
\arguments{
\item{IC}{a vector of IC values}
}
\value{
a vector of IC based weights
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

AICc.w <- weight.IC(IC=PGLSi$optim$AICc)

}
