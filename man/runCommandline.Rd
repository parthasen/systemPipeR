\name{runCommandline}
\alias{runCommandline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Execute SYSargs
}
\description{
Function to execute system parameters specified in \code{SYSargs} object
}
\usage{
runCommandline(args, runid = "01")
}
\arguments{
  \item{args}{
	object of class \code{SYSargs}
}
  \item{runid}{
	Run identifier used for log file to track system call commands. Default is \code{"01"}.
}
}
\details{
...
}
\value{
	Output files, their paths can be obtained with \code{outpaths()} from \code{SYSargs} container. 
	In addition, a \code{character vector} is returned containing the same paths.
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

\dontrun{
## Execute SYSargs on single machine
runCommandline(args=args)

## Execute SYSargs on multiple machines  
qsubargs <- getQsubargs(queue="batch", Nnodes="nodes=1", cores=cores(tophat), memory="mem=10gb", time="walltime=20:00:00")
qsubRun(args=args, qsubargs=qsubargs, Nqsubs=1, package="systemPipeR")
## Alignment stats
read_statsDF <- alignStats(fqpaths=tophatargs$infile1, bampaths=bampaths, fqgz=TRUE) 
read_statsDF <- cbind(read_statsDF[targets$FileName,], targets)
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
}
}
\keyword{ utilities }

