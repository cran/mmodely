\name{get.phylo.stats}
\alias{get.phylo.stats}
\title{Get tree statistics for a trait}
\description{
This function uses Pagel's lambda, Blombergs k, and Ancestral Character Estimation [ACE] to calculate statistics on a tree given a specified trait.
}
\usage{
get.phylo.stats(phylo, data, trait.clmn, gs.clmn='gn_sp', 
                ace.method='REML',ace.scaled=TRUE, ace.kappa=1)
}
\arguments{
\item{phylo}{PARAMDESCRIPTION}
\item{data}{PARAMDESCRIPTION}
\item{trait.clmn}{PARAMDESCRIPTION}
\item{gs.clmn}{PARAMDESCRIPTION}
\item{ace.method}{PARAMDESCRIPTION}
\item{ace.scaled}{PARAMDESCRIPTION}
\item{ace.kappa}{PARAMDESCRIPTION}
}
\value{
statistics on a particular trait within a tree (Pagel's lambda, Blomberg's  k, and the most ancestral ACE estimate) 
}
\examples{

data.path <- system.file("extdata","primate-example.data.csv", package="mmodely")
data <- read.csv(data.path, row.names=1)

data$gn_sp <- rownames(data)

tree.path <- system.file("extdata","primate-springer.2012.tre", package="mmodely")
phyl <- ape::read.tree(tree.path)[[5]]

get.phylo.stats(phylo=phyl, data=data, trait.clmn='OC', 
        gs.clmn='gn_sp', ace.method='REML',ace.scaled=TRUE, ace.kappa=1)

}
