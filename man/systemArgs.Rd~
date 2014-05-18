\name{systemArgs}
\alias{systemArgs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Command-line argument definitions
}
\description{
Return argurments of supported command-line software in \code{list}. 
}
\usage{
systemArgs(app = "tophat2", mymodules, mydir, myargs, myref, mygff, mytargets, myindir = "/data/", myoutdir = "/results/")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{app}{
	Name of command-line software provided as character vector of length one. Currently, supported tools are: Bowtie 2 and TopHat 2 
}
  \item{mymodules}{
	Version names of command-line software and other required tools to be loaded from a module system if available. One or many software versions can be assigened here to \code{mymodules} as character vector of any length, e.g. \code{c("bowtie2/2.1.0", "tophat/2.0.8b")}.
}
  \item{mydir}{
	Directory from where systemPipeR will be executed. Can be obtained with \code{getwd()}. Commonly, this directory will contain the subdirectories \code{data} and \code{results}.
}
  \item{myargs}{
	Named character vector containing the required arguments of command-line tool, e.g. \code{c(software="tophat", p="-p 4", i="-i 30", I="-I 3000")}.
}
  \item{myref}{
	Path to indexed reference genome as character vector of length one.
}
  \item{mygff}{
	Path to GTF or GFF3 genome annotation file. Only required for certain RNA-Seq aligners such as TopHat 2.
}
  \item{mytargets}{
	Name of targets file to use. For its format, see example below. 
}
  \item{myindir}{
	Name of directory containing inputs such as reference genome of GTF/GFF3 annotation files. Default is \code{data}. FASTQ file locations are usually specified in targets file.
}
  \item{myoutdir}{
	Name of directory where outputs should be written to, e.g. BAM files. Default is \code{result}.	
}
}
\details{
...
}
\value{
	\code{list}
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
\dontrun{
## Generate input targets file. Note: for \code{qsubRun()} the file targets_run.txt 
## needs to contain absolute paths to FASTQ files in "FileName' column.
library(systemPipeR)
targets <- read.delim("targets.txt", comment.char = "#")
write.table(targets, "targets_run.txt", row.names=FALSE, quote=FALSE, sep="\t")

## Build Bowtie 2 index
system("bowtie2-build ./data/mygenome.fa ./data/bowtie2index/mygenome")

## Run as single process without submitting to cluster, e.g. via qsub -I
mymodules <- c("bowtie2/2.1.0", "tophat/2.0.8b")
myargs <- c(software="tophat", p="-p 4", g="-g 1", segment_length="--segment-length 25", i="-i 30", I="-I 3000")
myref <- "./data/My_genome.fasta"
tophatargs <- systemArgs(app="tophat2", mymodules=mymodules, mydir=getwd(), myargs=myargs, myref=myref, mygff="My_specie.gff", mytargets="targets_run.txt")
bampaths <- runTophat(tophatargs=tophatargs, runid="01")

## Submit to compute nodes
qsubargs <- getQsubargs(queue="batch", Nnodes="nodes=4", cores=as.numeric(gsub("^.* ", "", tophatargs$args["p"])), memory="mem=10gb", time="walltime=20:00:00")
(joblist <- qsubRun(appfct="runTophat(appargs, runid)", appargs=tophatargs, qsubargs=qsubargs, Nqsubs=6, submitdir="/bigdata/tgirke/Projects/project_name/RNAseq/data", myfct="~/Projects/project_name/RNA-Seq/systemPipe.R"))

## Alignment stats
read_statsDF <- alignStats(fqpaths=tophatargs$infile1, bampaths=bampaths, fqgz=TRUE) 
read_statsDF <- cbind(read_statsDF[targets$FileName,], targets)
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")

}
}
\keyword{ utilities }
