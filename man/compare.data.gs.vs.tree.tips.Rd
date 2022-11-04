\name{compare.data.gs.vs.tree.tips}
\alias{compare.data.gs.vs.tree.tips}
\title{Find data being dropped by mismatches to the tree}
\description{
This function simply lists the rows of the data that are not getting matched to tips of the tree.
}
\usage{
compare.data.gs.vs.tree.tips(data, phylo, match.on=c('gn_sp','rownames')[1])
}
\arguments{
\item{data}{a data frame with genus species information as row names and a column named "gn_sp" }
\item{phylo}{a phylogenetic tree with labeled tip}
\item{match.on}{use a character string specifiying where the 'Genus_species' vector lies}
}
\value{
prints rows that are not matched ot the tree tips
}
\examples{

data.path <- system.file("extdata","primate-example.data.csv", package="mmodely")
data <- read.csv(data.path, row.names=1)
data$gn_sp <- rownames(data)

tree.path <- system.file("extdata","primate-springer.2012.tre", package="mmodely")
phyl <- ape::read.tree(tree.path)[[5]]


compare.data.gs.vs.tree.tips(data, phyl, match.on='rownames')

}
