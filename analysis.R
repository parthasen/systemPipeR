####################################
## Alignment with Bowtie2/Tophat2 ##
####################################
## Build Bowtie2 Index
system("bowtie2-build ./data/mygenome.fa ./data/bowtie2index/mygenome")

## Generate input targets file. Note: the name of input targets file is currently hard coded as: "targets_run.txt"!!!
targets <- read.delim("targets.txt", comment.char = "#")
write.table(targets[1,], "targets_run.txt", row.names=FALSE, quote=FALSE, sep="\t")

## Run as single process without submitting to cluster, e.g. in interactive session on owl or with qsub -I
source("~/Projects/project_name/RNA-Seq/systemPipe.R")
bampaths <- runTophat(tophatargs=tophatargs, runid="01")

## Submit to compute nodes
source("~/Projects/project_name/RNA-Seq/systemPipe.R")
(joblist <- qsubRun(appfct="runTophat(appargs, runid)", appargs=tophatargs, qsubargs=qsubargs, Nqsubs=6, submitdir="/bigdata/tgirke/Projects/project_name/RNA-Seq/results", myfct="~/Projects/project_name/RNA-Seq/systemPipe.R"))

## Alignment Stats
read_statsDF <- alignStats(fqpaths=tophatargs$infile1, bampaths=names(bampaths), fqgz=TRUE) 
read_statsDF <- cbind(read_statsDF[targets$FileName,], targets)
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
