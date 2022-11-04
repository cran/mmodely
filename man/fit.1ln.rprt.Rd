\name{fit.1ln.rprt}
\alias{fit.1ln.rprt}
\title{Report a model fit in a single line of text output}
\description{
This function takes a fit multivariate regression model as input and converts the normal tabular output into a single line using repeated "+"or"-" symbols for significance
}
\usage{
fit.1ln.rprt(fit, method=c('std.dev','p-value')[1], decimal.places=3, 
             name.char.len=6, print.inline=TRUE, rtrn.line=FALSE, R2AIC=TRUE,mn='')
}
\arguments{
\item{fit}{a fit model}
\item{method}{how to calculate the number of pluses or minuses before each coefficient name (default is standard deviations)}
\item{decimal.places}{the number of decimal places to use in reporting p-values}
\item{name.char.len}{the maximum length to use when truncating variable names}
\item{R2AIC}{boolean for also returning/printing AIC and R^2 values}
\item{print.inline}{should the outout string be printed to the terminal?}
\item{rtrn.line}{should the output string be returned as a characters string?}
\item{mn}{model number prefixed to printout if 'print.inline' is TRUE}
}
\value{
A character string of the form "++var1 +var5 var3 | -var2 --var4" indicating signifcance and direction of regression results
}
\examples{

path <- system.file("extdata","primate-example.data.csv", package="mmodely")
data <- read.csv(path, row.names=1)

model.fit <- lm('OC ~ mass.Kg + group.size + arboreal + leap.pct', data=data)

fit.1ln.rprt(fit=model.fit, decimal.places=3, name.char.len=6, print.inline=TRUE, rtrn.line=FALSE)

}
