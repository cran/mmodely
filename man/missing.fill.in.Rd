\name{missing.fill.in}
\alias{missing.fill.in}
\title{Fill in missing values in a dataframe with a secondary source}
\description{
This function uses the (non-missing) values from one column to fill in the missing values of another 
}
\usage{
missing.fill.in(x, var.from, var.to)
}
\arguments{
\item{x}{a dataframe or matrix}
\item{var.from}{secondary variable (of the same type and units) providing values to 'var.to'}
\item{var.to}{primary variable with missing values to fill in by 'var.from'}
}
\value{
a modified dataframe with fewer missing values in the 'var.to' column
}
\examples{

df <- data.frame(a=c(1,2,NA),b=c(1,NA,3),c=c(1,2,6))
missing.fill.in(df, 'c','a')

}
