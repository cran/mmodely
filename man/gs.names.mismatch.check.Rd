\name{gs.names.mismatch.check}
\alias{gs.names.mismatch.check}
\title{Check "Genus species" name formatting}
\description{
This convienience function checks to make sure that all of the elements the provided character vector adhere to the "Genus species" naming convention format. Default delimiters between genus and species names in the string are " ", "_", or "." 
}
\usage{
gs.names.mismatch.check(df, alias.table.path, gs.clmn='gn_sp')
}
\arguments{
\item{df}{a data frame with genus species information as row names and optionally in a column named "gn_sp" }
\item{alias.table.path}{a file system path (e.g. 'inst/extdata/primate.taxa.aliases.tab') to a lookup table with 'old.name' and 'new.name' as columns}
\item{gs.clmn}{the name of the column containing the 'Genus_species' vector}
}
\value{
None
}
\examples{

path <- system.file("extdata","primate-example.data.csv", package="mmodely")
gs.tab <- read.csv(path, row.names=1)
gs.tab$gn_sp <- rownames(gs.tab)

path.look <- system.file("extdata","primate.taxa.aliases.tab", package="mmodely")

gs.names.mismatch.check(gs.tab, alias.table.path=path.look, gs.clmn='gn_sp')

}
