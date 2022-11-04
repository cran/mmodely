\name{plot.confound.grid}
\alias{plot.confound.grid}
\title{Plot a grid of x y plots split by a confounder z}
\description{
Plot a grid of x y plots showing how a third confounding variable 'z' changes the slope
}
\usage{
\method{plot}{confound.grid}(x, Y='y', X='x', confounder='z', breaks=3,...)
}
\arguments{
\item{x}{a data frame}
\item{Y}{the name of the column with the dependent/outcome variable}
\item{X}{the name of the column with the predictor variable}
\item{confounder}{the name of the column with confounding variable}
\item{breaks}{number or vector of breaks to split the plots horizontally (across x)}
\item{...}{other arguments passed to 'plot'}
}
\value{
a confound grid plot
}
\examples{

path <- system.file("extdata","primate-example.data.csv", package="mmodely")
data <- read.csv(path, row.names=1)
data$col <- c('yellow','red')[data$nocturnal+1]

plot.confound.grid(x=data, Y='OC', X='leap.pct', confounder='mass.Kg')

}
