\name{plot.pgls.iters}
\alias{plot.pgls.iters}
\title{Plot the PGLS iterations}
\description{
A plot of AIC (and AICc) vs R^2 (and adjusted R^2) for all of the PGLS iterations
}
\usage{
\method{plot}{pgls.iters}(x, 
      bests=bestBy(x$optim, by=c('n','q','qXn','rwGsm')[1], best=c('AICc','R2.adj')[1], 
                         inverse=FALSE), ...)
}
\arguments{
\item{x}{a PGLSi[teration] object (a list of pgls model fits as well as optimization and tree parameter tables)}
\item{bests}{a table of the 'best' models to highlight in the plot based on some optimization criterion (e.g. R2)}
\item{...}{other parameters passed to 'plot'}
}
\value{
a plot of all of PGLS iterations
}
\examples{

data.path <- system.file("extdata","primate-example.data.csv", package="mmodely")
data <- read.csv(data.path, row.names=1)
pvs <- names(data[3:5])
data$gn_sp <- rownames(data)

tree.path <- system.file("extdata","primate-springer.2012.tre", package="mmodely")
phyl <- ape::read.tree(tree.path)[[5]]

mods <- get.model.combos(predictor.vars=pvs, outcome.var='OC', min.q=2)

PGLSi <- pgls.iter(models=mods, phylo=phyl, df=data, k=1,l=1,d=1) 


# sprinkle in some missing data so as to make model selection more interesting
for(pv in pvs){ data[sample(x=1:nrow(data),size=2),pv] <- NA} 

PGLSi <- pgls.iter(models=mods, phylo=phyl, df=data, k=1,l=1,d=1) 

# find the lowest AIC within each q by n sized sub-dataset
plot.pgls.iters(x=PGLSi)

}
