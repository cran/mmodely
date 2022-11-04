\name{comp.data}
\alias{comp.data}
\title{Comparative Data}
\description{
This is a shortcut function that wraps around "comparative.data" for use in the PGLS function.
}
\usage{
comp.data(phylo,df,gn_sp='gn_sp')
}
\arguments{
\item{phylo}{a tre file of the format phylo}
\item{df}{a data.frame with row names matching number and tip labels of 'tree'}
\item{gn_sp}{the column name (e.g. "gn_sp") that indicates how to match df with tree}
}
\value{
a "comparative data" table for use in PGLS modeling
}
\examples{

data.path <- system.file("extdata","primate-example.data.csv", package="mmodely")
data <- read.csv(data.path, row.names=1)
pvs <- names(data[3:5])
data$gn_sp <- rownames(data)

tree.path <- system.file("extdata","primate-springer.2012.tre", package="mmodely")
phyl <- ape::read.tree(tree.path)[[5]]

comp <- comp.data(phylo=phyl, df=data)


}
