\name{weight.AIC}
\alias{weight.AIC}
\title{Get AIC weights }
\description{
An implementation of AIC weighting that first calulates the difference in AIC values (here an AICc) by subtracting all values from the lowest AIC value. Second, the changes are expoentiated divided by a sum of the same and exponentiated yet again.
}
\usage{
weight.AIC(AIC)
}
\arguments{
\item{AIC}{a vector of [corrected] AIC values [AICc]}
}
\value{
a vector of AIC based weights
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

AICc.w <- weight.AIC(AIC=PGLSi$optim$AICc)

}
