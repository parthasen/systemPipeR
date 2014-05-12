##########################
## Environment Settings ##
##########################

## Load packages and functions
library(BSgenome); library(Rsamtools); library(rtracklayer); library(GenomicFeatures); library(Gviz); library(parallel); library(BiocParallel)
source("systemPipe.R")

## Generate input targets file. Note: for 'qsubRun()' the file targets_run.txt needs to contain absolute paths to FASTQ files in the "FileName' column.
targets <- read.delim("targets.txt", comment.char = "#")
write.table(targets, "targets_run.txt", row.names=FALSE, quote=FALSE, sep="\t")

############################
## Alignment with Tophat2 ##
############################
## Build Bowtie2 Index
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

## Alignment Stats
read_statsDF <- alignStats(fqpaths=tophatargs$infile1, bampaths=bampaths, fqgz=TRUE) 
read_statsDF <- cbind(read_statsDF[targets$FileName,], targets)
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")

##################################################################
## Alignment with Bowtie2 (here for miRNA profiling experiment) ##
##################################################################
## Run as single process without submitting to cluster, e.g. via qsub -I
mymodules <- c("bowtie2/2.1.0")
myargs <- c(software="bowtie2", p="-p 4", k="-k 50", other="--non-deterministic")
myref <- "./data/My_genome.fasta"
bowtieargs <- systemArgs(app="bowtie2", mymodules=mymodules, mydir=getwd(), myargs=myargs, myref=myref, mytargets="targets_run.txt")
bampaths <- runBowtie(bowtieargs=bowtieargs, runid="01")

## Submit to compute nodes
qsubargs <- getQsubargs(queue="batch", Nnodes="nodes=4", cores=as.numeric(gsub("^.* ", "", bowtieargs$args["p"])), memory="mem=10gb", time="walltime=20:00:00")
(joblist <- qsubRun(appfct="runBowtie(appargs, runid)", appargs=tophatargs, qsubargs=qsubargs, Nqsubs=6, submitdir="/bigdata/tgirke/Projects/project_name/RNAseq/data", myfct="~/Projects/project_name/RNA-Seq/systemPipe.R"))

##################################################
## Read counting for mRNA profiling experiments ##
##################################################
## Create txdb (do only once)
txdb <- makeTranscriptDbFromGFF(file="data/mygenome.gtf", format="gtf", dataSource="ENSEMBL", species="My_species")
saveDb(txdb, file="./data/My_species.sqlite")
## Read counting with summarizeOverlaps in parallel mode with multiple cores
library(BiocParallel)
txdb <- loadDb("./data/My_species.sqlite")
eByg <- exonsBy(txdb, by="gene")
bams <- names(bampaths); names(bams) <- targets$SampleName
bfl <- BamFileList(bams, yieldSize=50000, index=character())
multicoreParam <- MulticoreParam(workers=4); register(multicoreParam); registered()
countDFeByg <- bplapply(bfl, function(x) summarizeOverlaps(gff, x, mode="Union", ignore.strand=TRUE, inter.feature=TRUE, singleEnd=TRUE)) # Note: for strand-specific RNA-Seq set 'ignore.strand=FALSE' and for PE data set 'singleEnd=FALSE'
countDFeByg <- as.data.frame(row.names=names(rowData(countDFeByg[[1]])), sapply(seq(along=countDFeByg), function(x) assays(countDFeByg[[x]])$counts))
colnames(countDFeByg) <- names(bfl)
rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, gffsub=eByg))
write.table(assays(countDFeByg)$counts, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")

###################################################
## Read counting for miRNA profiling experiments ##
###################################################
## Download miRNA genes from miRBase
system("wget ftp://mirbase.org/pub/mirbase/19/genomes/My_species.gff3 -P ./data/")
gff <- import.gff("./data/My_species.gff3", asRangedData=FALSE)
gff <- split(gff, elementMetadata(gff)$ID)
bams <- names(bampaths); names(bams) <- targets$SampleName
bfl <- BamFileList(bams, yieldSize=50000, index=character())
countDFmiR <- summarizeOverlaps(gff, bfl, mode="Union", ignore.strand=FALSE, inter.feature=FALSE) # Note: inter.feature=FALSE important since pre and mature miRNA ranges overlap
rpkmDFmiR <- apply(countDFmiR, 2, function(x) returnRPKM(counts=x, gffsub=gff))
write.table(assays(countDFmiR)$counts, "results/countDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")
write.table(rpkmDFmiR, "results/rpkmDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")

#####################################
## Correlation analysis of samples ##
#####################################
library(ape)
rpkmDFeByg <- read.table("./results/rpkmDFeByg.xls", check.names=FALSE)
rpkmDFeByg <- rpkmDFeByg[rowMeans(rpkmDFeByg) > 50,]
d <- cor(rpkmDFeByg, method="spearman")
hc <- hclust(as.dist(1-d))
plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)



