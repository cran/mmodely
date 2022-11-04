\name{drop.na.data}
\alias{drop.na.data}
\title{Drop any rows with NA values}
\description{
This function takes a dataframe as input and removes any rows that have NA as values.
}
\usage{
drop.na.data(df, vars=names(df))
}
\arguments{
\item{df}{a dataframe}
\item{vars}{sub set of variable (column) names to use in searching for missing values}
}
\value{
A subset of 'df' that only has non-missing values in the columns specified by 'vars'
}
\examples{

path <- system.file("extdata","primate-example.data.csv", package="mmodely")
data <- read.csv(path, row.names=1)

df.nona <- drop.na.data(data, vars=names(df))


}
