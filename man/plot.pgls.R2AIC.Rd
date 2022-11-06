\name{plot.pgls.R2AIC}
\alias{plot.pgls.R2AIC}
\title{Plot (R2 vs AIC) results of a collection of fit PGLS models}
\description{
Plots a single panel of R^2 versus AIC, using versions of your choosing.
}
\usage{
\method{plot}{pgls.R2AIC}(x, 
            bests=bestBy(x, by=c('n','q','qXn','rwGsm')[4], best=c('AICc','R2.adj')[1], 
inverse=c(FALSE,TRUE)[1]), model.as.title='', bcl=rgb(1,1,1,maxColorValue=3,alpha=1),...)
}
\arguments{
\item{x}{a PGLSi[teration]$optim [optimization] table }
\item{bests}{a list of the best PGLS models grouped by variable count and sorted by some metric (e.g. adjusted R2)}
\item{model.as.title}{uses model.1ln.report to create a short character string of the "best" model results as a title}
\item{bcl}{background color of plot point}
\item{...}{other parameters passed to 'plot'}
}
\value{
a plot of R2 versus AIC of many PGLS models
}
\examples{

data.path <- system.file("extdata","primate-example.data.csv", package="mmodely")
data <- read.csv(data.path, row.names=1)
pvs <- names(data[3:6])
data$gn_sp <- rownames(data)

tree.path <- system.file("extdata","primate-springer.2012.tre", package="mmodely")
phyl <- ape::read.tree(tree.path)[[5]]

mods <- get.model.combos(predictor.vars=pvs, outcome.var='OC', min.q=2)

# sprinkle in some missing data so as to make model selection more interesting
for(pv in pvs){ data[sample(x=1:nrow(data),size=2),pv] <- NA} 

PGLSi <- pgls.iter(models=mods, phylo=phyl, df=data, k=1,l=1,d=1) 

plot.pgls.R2AIC(PGLSi$optim)   # find the lowest AIC within each q by n sized sub-dataset

}
