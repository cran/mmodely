\name{gs.rename}
\alias{gs.rename}
\title{Rename the Genus species information in a data frame}
\description{
This function takes a data frame (with a genus species column) and proceeds to use an external look-up table to update the names if they've been changed
}
\usage{
gs.rename(df, alias.table.path, retro=FALSE, update.gn_sp=FALSE)
}
\arguments{
\item{df}{a data frame with genus species information as row names and optionally in a column named "gn_sp" }
\item{alias.table.path}{a file system path (e.g. 'inst/extdata/primate.taxa.aliases.tab') to a lookup table with 'old.name' and 'new.name' as columns}
\item{retro}{a boolean (T/F) parameter specifying if the renaming should go from new to old instead of the default of old to new}
\item{update.gn_sp}{a boolean parameter specifying if the 'gn_sp' column should also be updated with 'new.name's}
}
\value{
the original data frame with (potentially) updated row names and updated gn_sp column values
}
\examples{

path.data <- system.file("extdata","primate-example.data.csv", package="mmodely")
data <- read.csv(path.data, row.names=1)

path.look <- system.file("extdata","primate.taxa.aliases.tab", package="mmodely")

data.renamed <- gs.rename(df=data, alias.table.path=path.look, retro=FALSE, update.gn_sp=FALSE)

}
