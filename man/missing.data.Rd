\name{missing.data}
\alias{missing.data}
\title{Report missing values in a dataframe}
\description{
This funciton reports column and rowwise missing data. It can also list the rownames for missing columns or the column names for missing rows.
}
\usage{
missing.data(x, cols=NULL, rows=NULL)
}
\arguments{
\item{x}{a dataframe}
\item{cols}{print the specific rows corresponding to missing values in this column}
\item{rows}{print the specific cols corresponding to missing values in this rowname}
}
\value{
a report on column versus row wise missing data
}
\examples{

path <- system.file("extdata","primate-example.data.csv", package="mmodely")
data <- read.csv(path, row.names=1)

missing.data(data)


}
