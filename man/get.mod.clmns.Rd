\name{get.mod.clmns}
\alias{get.mod.clmns}
\title{Get model columns}
\description{
Get the variable names from a model string by splitting on "+" and '~' using both 'get.mod.outcome' and 'get.mod.vars'. The results are passed to the comp.data function for eventual use in PGLS modeling. 'gn_sp' is included as it is typically required to link tree tips to rows of the comparative data. 
}
\usage{
get.mod.clmns(model, gs.clmn='gn_sp')
}
\arguments{
\item{model}{a model string of the form "y ~ x1 + x2 ..."}
\item{gs.clmn}{the column header for the vector of "Genus_species" names, to link a tree tips to rows}
}
\value{
a vector of characters enummerating the columns to retain in PGLS modeling (input to df param in the 'comp.data' function)
}
\examples{

model.columns <- get.mod.clmns(model=formula('y ~ x1 + x2'))

}
