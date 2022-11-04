

\name{modsel.distro.dots}
\alias{modsel.distro.dots}
\title{Coeficients distribution plot of models selected from each subset}
\description{
Plot the raw distribution of points corresponding to the coefficients harvested from the best model of each subset of the dataset.
}
\usage{
modsel.distro.dots(PC, jit.f=1, R2x=3, zeroline=TRUE, add=FALSE, pd=0, 
                  pvs=names(PC$coefs), pvlabs=NULL, xlim=range(unlist(PC$coefs)),...)
}
\arguments{
\item{PC}{a list of vectors of pooled coefficients (or scores) harvested from the 'best' selected modeling runs (out put from 'get.pgls.coefs') }
\item{jit.f}{factor for random jittering (see 'jitter()'}
\item{R2x}{the character expansion factor according to R^2 value}
\item{zeroline}{should we add an abline at x=0?}
\item{add}{should we add to the existing plot?}
\item{pd}{'position dodge' moves all y axis plotting positions up or down by this provided value (useful for adding multiple distributions for the same param)}
\item{pvs}{the predictor variable vector for ordering the y-axis labels}
\item{pvlabs}{the predictor variable labels for labeling the plot (defaults to pvs)}
\item{xlim}{x axis plot limits}
\item{...}{other parameters passed on to plot}
}
\value{
a 'distro dot plot' of coefficent distributions
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

a.PGLS <- get.best.model(PGLSi, by=c('R2.adj','AICc')[1])

coefs.objs <- get.pgls.coefs(PGLSi$fits, est='Estimate')

modsel.distro.dots(coefs.objs)

}
