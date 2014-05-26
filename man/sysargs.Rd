\name{sysargs}
\alias{modules}
\alias{modules-methods}
\alias{modules,SYSargs-method}
\alias{software}
\alias{software-methods}
\alias{software,SYSargs-method}
\alias{cores}
\alias{cores-methods}
\alias{cores,SYSargs-method}
\alias{other}
\alias{other-methods}
\alias{other,SYSargs-method}
\alias{reference}
\alias{reference-methods}
\alias{reference,SYSargs-method}
\alias{infile1}
\alias{infile1-methods}
\alias{infile1,SYSargs-method}
\alias{infile2}
\alias{infile2-methods}
\alias{infile2,SYSargs-method}
\alias{outfile1}
\alias{outfile1-methods}
\alias{outfile1,SYSargs-method}
\alias{SampleName}
\alias{SampleName-methods}
\alias{SampleName,SYSargs-method}
\alias{sysargs-methods}
\alias{sysargs,SYSargs-method}
\alias{outpaths}
\alias{outpaths-methods}
\alias{outpaths,SYSargs-method}
\docType{data}
\title{ SYSargs accessor methods}

\description{
	Methods to access information from \code{SYSargs} object.
}
\usage{
	sysargs(x)
}
\arguments{
  \item{x}{
	object of class \code{SYSargs}
}
}
\details{
...
}
\value{
	various outputs
}
\references{
...
}
\author{
Thomas Girke
}
\seealso{
...
}
\examples{
## Construct SYSargs object from param and targets files 
param <- paste0(system.file("extdata", package="systemPipeR"), "/tophat.param")
targets <- paste0(system.file("extdata", package="systemPipeR"), "/targets.txt")
args <- systemArgs(sysma=param, mytargets=targets)
args
names(args); modules(args); cores(args); outpaths(args); sysargs(args)
}
\keyword{ utilities }

