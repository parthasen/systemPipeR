\name{getQsubargs}
\alias{getQsubargs}
\title{
Arguments for qsub 
}
\description{
Defines arguments to submit runX job(s) to queuing system (e.g. Torque) via \code{qsub}.
}
\usage{
getQsubargs(software = "qsub", queue = "batch", Nnodes = "nodes=1", cores = as.numeric(gsub("^.* ", "", tophatargs$args["p"])), memory = "mem=10gb", time = "walltime=20:00:00")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{software}{
	Software to use for submission to queuing system. Default is \code{qsub}.
}
  \item{queue}{
	Name of queue to to use. Default is \code{batch}.
}
  \item{Nnodes}{
	Number of compute nodes to use for processing. Default is \code{nodes=1}.
}
  \item{cores}{
	Number of CPU cores to use per compute node. Default will use what is provided by under \code{-p} in \code{myargs} of \code{systemArgs()} output.
}
  \item{memory}{
	Amount of RAM to reserve per node.
}
  \item{time}{
	Walltime limit each job is allowed to run per node.
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
