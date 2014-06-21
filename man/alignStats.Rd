\name{alignStats}
\alias{alignStats}
\title{
Alignment statistics 
}
\description{
Generate data frame containg important read alignment statistics such as the
total number of reads in the FASTQ files, the number of total alignments, as
well as the number of primary alignments in the corresponding BAM files. 
}
\usage{
alignStats(args, fqgz = TRUE)
}
\arguments{
  \item{args}{
	Object of class \code{SYSargs}. 
}
  \item{fqgz}{
	Defines whether FASTQ files are compressed (\code{TRUE} or \code{FALSE}).
}
}
\details{
...
}
\value{
	\code{data.frame} with alignment statistics.
}
\references{
...
}
\author{
Thomas Girke
}
\note{
...
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
read_statsDF <- alignStats(args, fqgz=TRUE) 
read_statsDF <- cbind(read_statsDF[targets$FileName,], targets)
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
}
}
\keyword{ utilities }
