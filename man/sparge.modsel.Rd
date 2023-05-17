\name{sparge.modsel}
\alias{sparge.modsel}
\title{Coeficients distribution [sparge] plot of models selected from each subset}
\description{
Plot the raw distribution of points corresponding to the coefficients harvested from the best model of each subset of the dataset.
}
\usage{
sparge.modsel(PC, jit.f=1, R2x=3, nx=2, n.max=max(unlist(PC$n)), zeroline=TRUE,
                  add=FALSE, pd=0, pvs=names(PC$coefs), pvlabs=NULL, 
                  xlim=range(unlist(PC$coefs)),
                  MA = NULL, ap=8, ac = 1, ax = nx, ...)
}
\arguments{
\item{PC}{a list of vectors of pooled coefficients (or scores) harvested from the 'best' selected modeling runs (out put from 'get.pgls.coefs') }
\item{jit.f}{factor for random jittering (see 'jitter()'}
\item{R2x}{the line width expansion factor according to R^2 value}
\item{nx}{the point size expansion factor according to sample size of model} 
\item{n.max}{the maximum sample size used in all models} 
\item{zeroline}{should we add an abline at x=0?}
\item{add}{should we add to the existing plot?}
\item{pd}{'position dodge' moves all y axis plotting positions up or down by this provided value (useful for adding multiple distributions for the same param)}
\item{pvs}{the predictor variable vector for ordering the y-axis labels}
\item{pvlabs}{the predictor variable labels for labeling the plot (defaults to pvs)}
\item{xlim}{x axis plot limits}
\item{MA}{matrix of model averages (defaults to NULL)}
\item{ap}{coded numeric point character symbol used for model averaged parameter position}
\item{ac}{color symbol used for model averaged parameters plot character}
\item{ax}{expansion factor to expant model average parameter plot character (defaults to nx)}
\item{...}{other parameters passed on to plot}
}
\value{
a 'sparge' [sprinkle/smear] plot of coefficent distributions
}
\seealso{
See also 'boxplot' and 'stripchart' in package 'graphics' as well as 'violin', 'bean', 'ridgelines', and 'raincloud' plots.
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

coefs.objs <- get.pgls.coefs(PGLSi$fits, est='Estimate')

sparge.modsel(coefs.objs)

}
