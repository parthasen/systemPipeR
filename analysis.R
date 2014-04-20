####################################
## Alignment with Bowtie2/Tophat2 ##
####################################
## Build Bowtie2 Index
system("bowtie2-build ./data/mygenome.fa ./data/bowtie2index/mygenome")

## Generate input targets file. Note: for 'qsubRun()' the file targets_run.txt needs to contain absolute paths to FASTQ files in the "FileName' column.
targets <- read.delim("targets.txt", comment.char = "#")
write.table(targets[1,], "targets_run.txt", row.names=FALSE, quote=FALSE, sep="\t")

## Run as single process without submitting to cluster, e.g. in interactive session on owl or with qsub -I
source("systemPipe.R")
tophatargs <- systemArgs(mydir=getwd(), myref="TAIR10_chr_all.fas", mygff="TAIR10_GFF3_genes.gff", mytargets="targets_run.txt")
bampaths <- runTophat(tophatargs=tophatargs, runid="01")

## Submit to compute nodes
qsubargs <- getQsubargs(queue="batch", Nnodes="nodes=4", cores=as.numeric(gsub("^.* ", "", tophatargs$args["p"])), memory="mem=10gb", time="walltime=20:00:00")
(joblist <- qsubRun(appfct="runTophat(appargs, runid)", appargs=tophatargs, qsubargs=qsubargs, Nqsubs=6, submitdir="/bigdata/tgirke/Projects/project_name/RNAseq/data", myfct="~/Projects/project_name/RNA-Seq/systemPipe.R"))

## Alignment Stats
read_statsDF <- alignStats(fqpaths=tophatargs$infile1, bampaths=names(bampaths), fqgz=TRUE) 
read_statsDF <- cbind(read_statsDF[targets$FileName,], targets)
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
