\name{cept}
\alias{cept}
\title{Include all variables except ...}
\description{
This function takes a dataframe, list, or a named vector of variable (column) names to subset
}
\usage{
cept(x,except='gn_sp')
}
\arguments{
\item{x}{a dataframe, list, or named vector}
\item{except}{a vector of the names of the items in x to exclude}
}
\value{
the subset of x without those 'except' items specified
}
\examples{

data.path <- system.file("extdata","primate-example.data.csv", package="mmodely")
data <- read.csv(data.path, row.names=1)

df.except.gnsp <- cept(x=data,except='gn_sp')

}
