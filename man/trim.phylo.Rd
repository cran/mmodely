\name{trim.phylo}
\alias{trim.phylo}
\title{Trim a phylogenetic tree using Genus species names}
\description{
Read in a vector of genus species names and a tree and drop the tips in the tree that match the vector of names.
}
\usage{
trim.phylo(phylo, gs.vect)
}
\arguments{
\item{phylo}{a phylogenetic tree}
\item{gs.vect}{a vector of character strings in the 'Genus_species' format}
}
\value{
a plot of a transformed phylogenetic tree
}
\examples{

data.path <- system.file("extdata","primate-example.data.csv", package="mmodely")
data <- read.csv(data.path, row.names=1)
data$gn_sp <- rownames(data)

tree.path <- system.file("extdata","primate-springer.2012.tre", package="mmodely")
phylo <- read.tree(tree.path)[[5]]

trim.phylo(phylo, gs.vect=data$gn_sp)

}
