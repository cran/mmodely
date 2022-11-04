\name{plot.xy.ab.p}
\alias{plot.xy.ab.p}
\title{An x/y scatterplot with a linear regression line and p-value}
\description{
This function performs a simple scatter plot but also superimposses a linear regression trend (abline) and optionally also the p-value of this line
}
\usage{
\method{plot}{xy.ab.p}(x, x.var, y.var, 
fit.line=TRUE, p.value=TRUE, slope=TRUE, p.col='red', plot.labels=TRUE, verbose=TRUE, ...)
}
\arguments{
\item{x}{a data frame}
\item{x.var}{the name of the x variable in df}
\item{y.var}{the name of the y variable in df}
\item{fit.line}{should a fit (ab) line be drawn?}
\item{p.value}{should the p-value be printed on the plot?}
\item{slope}{should the slope be printed on the plot?}
\item{p.col}{should the plot be labeled?}
\item{plot.labels}{should all of thie model fit information be printed out?}
\item{verbose}{should all other information be printed out too?}
\item{...}{other parameters passed to 'plot'}
}
\value{
An x/y scatterplot with regression line
}
\examples{

path <- system.file("extdata","primate-example.data.csv", package="mmodely")
data <- read.csv(path, row.names=1)

plot.xy.ab.p(x=data, x.var='OC', y.var='group.size', 
        fit.line=TRUE, p.value=TRUE, slope=TRUE, p.col='red', plot.labels=TRUE, verbose=TRUE)

}
