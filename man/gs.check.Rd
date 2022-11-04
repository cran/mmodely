\name{gs.check}
\alias{gs.check}
\title{Check "Genus species" name formatting}
\description{
This convienience function checks to make sure that all of the elements the provided character vector adhere to the "Genus species" naming convention format. Default delimiters between genus and species names in the string are " ", "_", or "." 
}
\usage{
gs.check(genus.species, sep='[ _\\\\.]')
}
\arguments{
\item{genus.species}{a vector of character strings specifiying the combination of Genus [and] species}
\item{sep}{a regular expression between genus and species}
}
\value{
None
}
\examples{

path <- system.file("extdata","primate-example.data.csv", package="mmodely")
gs.tab <- read.csv(path, row.names=1)
gs.tab$gn_sp <- rownames(gs.tab)

gs.check(genus.species=gs.tab$gn_sp, sep='[ _\\\\.]')

}
