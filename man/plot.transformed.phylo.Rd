\name{plot.transformed.phylo}
\alias{plot.transformed.phylo}
\title{Plot a transformed phylogenetic tree}
\description{
PGLS regression will use maximum likelihood to estimate tree parameters while also estimating regression parameters. Here we provide a utility function to visualize what this new tree would look like in two dimensions.
}
\usage{
\method{plot}{transformed.phylo}(x, delta=1,kappa=1,...)
}
\arguments{
\item{x}{a phylogenetic tree}
\item{delta}{an integer between 0 and 3}
\item{kappa}{an integer between 0 and 3}
\item{...}{other parameters passed to 'plot'}
}
\value{
a plot of a transformed phylogenetic tree
}
\examples{

tree.path <- system.file("extdata","primate-springer.2012.tre", package="mmodely")
phyl <- ape::read.tree(tree.path)[[5]]

plot.transformed.phylo(x=phyl, delta=2.3,kappa=2.1)

}
