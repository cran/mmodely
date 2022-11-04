\name{interpolate}
\alias{interpolate}
\title{Interpolate missing data in a data frame}
\description{
This function finds NA values and interpolates using averaging values of nearby genus and species  
}
\usage{
interpolate(df, taxa=c('genus','family'), clmns=1:length(df))
}
\arguments{
\item{df}{a data frame}
\item{taxa}{a vector of taxonomic ranks (corresonding to columns) to assist in guiding the interpolating}
\item{clmns}{the names of the columns to interpolate over }
}
\value{
a modified data frame without missing values in the columns specified
}
\examples{

path <- system.file("extdata","primate-example.data.csv", package="mmodely")
gs.tab <- read.csv(path, row.names=1)

clmns <- match(c('mass.Kg','DPL.km'),names(gs.tab))
df.2 <- interpolate(df=gs.tab, taxa='genus', clmns=clmns)

}
