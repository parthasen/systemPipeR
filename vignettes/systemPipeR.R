### R code from vignette source 'systemPipeR.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: systemPipeR.Rnw:41-43
###################################################
options(width=95)
unlink("test.db")


###################################################
### code chunk number 3: systemPipeR.Rnw:65-67 (eval = FALSE)
###################################################
## # $ R CMD build systemPipeR # Builds package
## install.packages("systemPipeR.1.0.0.tar.gz", repos=NULL, type="source") # Installs the package


###################################################
### code chunk number 4: systemPipeR.Rnw:72-75 (eval = FALSE)
###################################################
## library("systemPipeR") # Loads the package
## library(help="systemPipeR") # Lists all functions and classes 
## vignette("systemPipeR") # Opens this PDF manual from R


###################################################
### code chunk number 5: systemPipeR.Rnw:80-83
###################################################
library(systemPipeR)
targetspath <- paste0(system.file("extdata", package="systemPipeR"), "/targets.txt")
read.delim(targetspath, comment.char = "#")


###################################################
### code chunk number 6: systemPipeR.Rnw:87-88
###################################################
readComp(file=targetspath, format="vector", delim="-")


###################################################
### code chunk number 7: systemPipeR.Rnw:94-96 (eval = FALSE)
###################################################
## library(systemPipeR)
## library(BSgenome); library(Rsamtools); library(rtracklayer); library(GenomicFeatures); library(Gviz); library(parallel); library(BiocParallel)


###################################################
### code chunk number 8: systemPipeR.Rnw:100-102 (eval = FALSE)
###################################################
## targets <- read.delim("targets.txt", comment.char = "#")
## write.table(targets, "targets_run.txt", row.names=FALSE, quote=FALSE, sep="\t")


###################################################
### code chunk number 9: systemPipeR.Rnw:107-108 (eval = FALSE)
###################################################
## system("bowtie2-build ./data/mygenome.fa ./data/bowtie2index/mygenome")


###################################################
### code chunk number 10: systemPipeR.Rnw:112-117 (eval = FALSE)
###################################################
## mymodules <- c("bowtie2/2.1.0", "tophat/2.0.8b")
## myargs <- c(software="tophat", p="-p 4", g="-g 1", segment_length="--segment-length 25", i="-i 30", I="-I 3000")
## myref <- "./data/My_genome.fasta"
## tophatargs <- systemArgs(app="tophat2", mymodules=mymodules, mydir=getwd(), myargs=myargs, myref=myref, mygff="My_specie.gff", mytargets="targets_run.txt")
## bampaths <- runTophat(tophatargs=tophatargs, runid="01")


###################################################
### code chunk number 11: systemPipeR.Rnw:121-123 (eval = FALSE)
###################################################
## qsubargs <- getQsubargs(queue="batch", Nnodes="nodes=4", cores=as.numeric(gsub("^.* ", "", tophatargs$args["p"])), memory="mem=10gb", time="walltime=20:00:00")
## (joblist <- qsubRun(appfct="runTophat(appargs, runid)", appargs=tophatargs, qsubargs=qsubargs, Nqsubs=6, submitdir="/bigdata/tgirke/Projects/project_name/RNAseq/data", myfct="~/Projects/project_name/RNA-Seq/systemPipe.R"))


###################################################
### code chunk number 12: systemPipeR.Rnw:127-130 (eval = FALSE)
###################################################
## read_statsDF <- alignStats(fqpaths=tophatargs$infile1, bampaths=bampaths, fqgz=TRUE) 
## read_statsDF <- cbind(read_statsDF[targets$FileName,], targets)
## write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")


###################################################
### code chunk number 13: systemPipeR.Rnw:135-140 (eval = FALSE)
###################################################
## mymodules <- c("bowtie2/2.1.0")
## myargs <- c(software="bowtie2", p="-p 4", k="-k 50", other="--non-deterministic")
## myref <- "./data/My_genome.fasta"
## bowtieargs <- systemArgs(app="bowtie2", mymodules=mymodules, mydir=getwd(), myargs=myargs, myref=myref, mytargets="targets_run.txt")
## bampaths <- runBowtie(bowtieargs=bowtieargs, runid="01")


###################################################
### code chunk number 14: systemPipeR.Rnw:144-146 (eval = FALSE)
###################################################
## qsubargs <- getQsubargs(queue="batch", Nnodes="nodes=4", cores=as.numeric(gsub("^.* ", "", bowtieargs$args["p"])), memory="mem=10gb", time="walltime=20:00:00")
## (joblist <- qsubRun(appfct="runBowtie(appargs, runid)", appargs=tophatargs, qsubargs=qsubargs, Nqsubs=6, submitdir="/bigdata/tgirke/Projects/project_name/RNAseq/data", myfct="~/Projects/project_name/RNA-Seq/systemPipe.R"))


###################################################
### code chunk number 15: systemPipeR.Rnw:151-153 (eval = FALSE)
###################################################
## txdb <- makeTranscriptDbFromGFF(file="data/mygenome.gtf", format="gtf", dataSource="ENSEMBL", species="My_species")
## saveDb(txdb, file="./data/My_species.sqlite")


###################################################
### code chunk number 16: systemPipeR.Rnw:157-169 (eval = FALSE)
###################################################
## library(BiocParallel)
## txdb <- loadDb("./data/My_species.sqlite")
## eByg <- exonsBy(txdb, by="gene")
## bams <- names(bampaths); names(bams) <- targets$SampleName
## bfl <- BamFileList(bams, yieldSize=50000, index=character())
## multicoreParam <- MulticoreParam(workers=4); register(multicoreParam); registered()
## counteByg <- bplapply(bfl, function(x) summarizeOverlaps(gff, x, mode="Union", ignore.strand=TRUE, inter.feature=TRUE, singleEnd=TRUE)) # Note: for strand-specific RNA-Seq set 'ignore.strand=FALSE' and for PE data set 'singleEnd=FALSE'
## countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
## rownames(countDFeByg) <- names(rowData(counteByg[[1]])); colnames(countDFeByg) <- names(bfl)
## rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, gffsub=eByg))
## write.table(assays(countDFeByg)$counts, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
## write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")


###################################################
### code chunk number 17: systemPipeR.Rnw:174-183 (eval = FALSE)
###################################################
## system("wget ftp://mirbase.org/pub/mirbase/19/genomes/My_species.gff3 -P ./data/")
## gff <- import.gff("./data/My_species.gff3", asRangedData=FALSE)
## gff <- split(gff, elementMetadata(gff)$ID)
## bams <- names(bampaths); names(bams) <- targets$SampleName
## bfl <- BamFileList(bams, yieldSize=50000, index=character())
## countDFmiR <- summarizeOverlaps(gff, bfl, mode="Union", ignore.strand=FALSE, inter.feature=FALSE) # Note: inter.feature=FALSE important since pre and mature miRNA ranges overlap
## rpkmDFmiR <- apply(countDFmiR, 2, function(x) returnRPKM(counts=x, gffsub=gff))
## write.table(assays(countDFmiR)$counts, "results/countDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")
## write.table(rpkmDFmiR, "results/rpkmDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")


###################################################
### code chunk number 18: systemPipeR.Rnw:187-193 (eval = FALSE)
###################################################
## library(ape)
## rpkmDFeByg <- read.table("./results/rpkmDFeByg.xls", check.names=FALSE)
## rpkmDFeByg <- rpkmDFeByg[rowMeans(rpkmDFeByg) > 50,]
## d <- cor(rpkmDFeByg, method="spearman")
## hc <- hclust(as.dist(1-d))
## plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)


###################################################
### code chunk number 19: systemPipeR.Rnw:197-198 (eval = FALSE)
###################################################
## cmp <- readComp(myfile="targets.txt", format="vector", delim="-")


###################################################
### code chunk number 20: sessionInfo
###################################################
toLatex(sessionInfo())


