\name{correct.AIC}
\alias{correct.AIC}
\title{Correct AIC}
\description{
Calculate a corrected Akaiki Information Criterion 
}
\usage{
correct.AIC(AIC, K,n)
}
\arguments{
\item{AIC}{a vector of AIC values}
\item{K}{number of parameters}
\item{n}{number of data}
}
\value{
corrected AIC values
}
\examples{

correct.AIC(AIC=100,K=10,n=100)

}
