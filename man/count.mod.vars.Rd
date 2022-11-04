\name{count.mod.vars}
\alias{count.mod.vars}
\title{Count the predictor variables in a model}
\description{
This function takes a model string and counts the number of predictor variables.
}
\usage{
count.mod.vars(model)
}
\arguments{
\item{model}{model specified as a string in the form "y ~ x1 + x2 ..."}
}
\value{
an integer specifying the count of predictor variables
}
\examples{

count <- count.mod.vars(model=formula('y ~ x1 + x2'))
if(count == 2) { print('sane'); }else{ print('insane')} 

}

