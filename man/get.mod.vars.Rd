\name{get.mod.vars}
\alias{get.mod.vars}
\title{Get model variable names}
\description{
Split the predictor string of a model formula into it's constituent character strings.
}
\usage{
get.mod.vars(model)
}
\arguments{
\item{model}{a character string of a formula of the form 'y ~ x1 + x2'}
}
\value{
a vector of character strings of variable names (e.g. corresponding to column names for comp.data input)
}
\examples{

model.variables <- get.mod.vars(model='y ~ x1 + x2')

}
