\name{get.model.combos}
\alias{get.model.combos}
\title{All combinations of predictor variables}
\description{
Enumerate all combinations of predictor variables in a multivariate regression model.
}
\usage{
get.model.combos(outcome.var, predictor.vars, min.q=2)
}
\arguments{
\item{predictor.vars}{predictor variables names (a vector of character strings)}
\item{outcome.var}{outcome variable name (character string)}
\item{min.q}{minimum number of predictor variables to include in the mode (default is 2) }
}
\value{
a vector of models as character strings of the form "y ~ x1 + x2 ..."
}
\examples{

path <- system.file("extdata","primate-example.data.csv", package="mmodely")
data <- read.csv(path, row.names=1)

get.model.combos(outcome.var='OC', predictor.vars=names(data),  min.q=2)

}
