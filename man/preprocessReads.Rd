\name{preprocessReads}
\alias{preprocessReads}
\title{
Run custom read preprocessing functions
}
\description{
Function to run custom read preprocessing functions on FASTQ files specified in
the \code{infile1} slot of \code{SYSargs} objects. The names of the corresponding 
output FASTQ files are specified in the \code{outpaths} slot of the same \code{SYSargs} 
object. The function uses the \code{FastqStreamer} function from the \code{ShortRead} 
package to stream through large files in a memory-efficient manner.
}
\usage{
preprocessReads(args, Fct, batchsize = 1e+05, overwrite = TRUE, ...)
}
\arguments{
  \item{args}{
Object of class \code{SYSargs}
}
  \item{Fct}{
\code{character} string of custom read preprocessing function call where both the 
input and output needs to be an object of class \code{ShortReadQ}. The name of the 
input \code{ShortReadQ} object needs to be \code{fq}.
}
  \item{batchsize}{
Number of reads to process in each iteration by the internally used \code{FastqStreamer}
function.
}
  \item{overwrite}{
If \code{TRUE} existing file will be overwritten.
}
  \item{\dots}{
To pass on additional arguments to the internally used \code{writeFastq} function.
}
}
\value{
Writes to files in FASTQ format. Their names are specified by \code{outpaths(args)}.
}
\author{
Thomas Girke
}
\seealso{
FastqStreamer 
}
\examples{
param <- system.file("extdata", "trim.param", package="systemPipeR")
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
args <- systemArgs(sysma=param, mytargets=targets)
\dontrun{
preprocessReads(args=args, Fct="trimLRPatterns(Rpattern='GCCCGGGTAA', subject=fq)", batchsize=100000, overwrite=TRUE, compress=TRUE)
}
}
\keyword{ utilities }
