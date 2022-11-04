\name{pgls.report}
\alias{pgls.report}
\title{Report PGLS results as a table}
\description{
Output a spreadsheet ready tabular summary of a fit PGLS model
}
\usage{
pgls.report(cd, f=formula('y~x'), l=1,k=1,d=1, 
            bounds=list(lambda=c(.2,1),kappa=c(.2,2.7),delta=c(.2,2.7)), 
            anova=FALSE, mod.no='NA', out='pgls.output-temp',QC.plot=FALSE)
}
\arguments{
\item{cd}{a comparative data object, here created by 'comp.data'}
\item{f}{the model formula (as acharacter string)}
\item{k}{the fixed or 'ML' value for kappa}
\item{l}{the fixed or 'ML' value for lambda}
\item{d}{the fixed or 'ML' value for delta}
\item{bounds}{a list of vectors of upper and lower bounds for kappa, lambda, and delta}
\item{anova}{should an anova be run on the fit model and output to the terminal?}
\item{mod.no}{the model number (can be the order that models were run)}
\item{out}{the base filename to be printed out}
\item{QC.plot}{should a quality control plot be output to screen?}
}
\value{
A summary results of a fit PGLS model with ANOVA and tabular spreadsheet ready csv filesystem output.
}
\examples{


data.path <- system.file("extdata","primate-example.data.csv", package="mmodely")
data <- read.csv(data.path, row.names=1)
data$gn_sp <- rownames(data)

tree.path <- system.file("extdata","primate-springer.2012.tre", package="mmodely")
#5. RAxML phylogram based on the 61199 bp concatenation of 69 nuclear and ten mitochondrial genes. 
phyl <- ape::read.tree(tree.path)[[5]]

phyl <- trim.phylo(phylo=phyl, gs.vect=data$gn_sp)

comp <- comp.data(phylo=phyl, df=data)

pgls.report(comp, f=formula('OC~mass.Kg + DPL.km'), l=1,k=1,d=1,
             anova=FALSE, mod.no='555', out='', QC.plot=TRUE)

}
