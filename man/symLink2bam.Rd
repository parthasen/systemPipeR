\name{symLink2bam}
\alias{symLink2bam}
\title{
Symbolic links for IGV
}
\description{
Function for creating symbolic links to view BAM files in a genome browser such as IGV.
}
\usage{
symLink2bam(sysargs, command="ln -s", htmldir, ext = c(".bam", ".bai"), urlbase, urlfile)
}
\arguments{
  \item{sysargs}{
	Object of class \code{SYSargs}
}
  \item{command}{
	Shell command, defaults to "ln -s"
}
  \item{htmldir}{
	Path to HTML directory with http access.
}
  \item{ext}{
	File name extensions to use for BAM and index files.
}
  \item{urlbase}{
	The base URL structure to use in URL file.
}
  \item{urlfile}{
	Name and path of URL file.
}
}
\details{
...
}
\value{
	symbolic links and url file
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

\dontrun{
## Create sym links and URL file for IGV
symLink2bam(sysargs=args, command="ln -s", htmldir=c("~/.html/", "somedir/"), ext=c(".bam", ".bai"), urlbase="http://myserver.edu/~username/", urlfile="IGVurl.txt")
}
}
\keyword{ utilities }

