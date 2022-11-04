\name{get.mod.outcome}
\alias{get.mod.outcome}
\title{Get the outcome variable from a model string}
\description{
Get the outcome variable from the front of a model formula string. Used as part of 'get.mod.clmns' to be passed to 'comp.data'
}
\usage{
get.mod.outcome(model)
}
\arguments{
\item{model}{a character string of a formula of the form 'y ~ x1 + x2 ...'}
}
\value{
a character string specifying the outcome variable
}
\examples{

model.columns <- get.mod.clmns(model=formula('y ~ x1 + x2'))

}
