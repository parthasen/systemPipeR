### R code from vignette source 'systemPipeR.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: systemPipeR.Rnw:42-44
###################################################
options(width=95)
unlink("test.db")


###################################################
### code chunk number 3: systemPipeR.Rnw:66-68 (eval = FALSE)
###################################################
## system("R CMD build systemPipeR") # Builds package
## install.packages("systemPipeR.X.X.X.tar.gz", repos=NULL, type="source") # Installs the package


###################################################
### code chunk number 4: systemPipeR.Rnw:73-76 (eval = FALSE)
###################################################
## library("systemPipeR") # Loads the package
## library(help="systemPipeR") # Lists all functions and classes 
## vignette("systemPipeR") # Opens this PDF manual from R


###################################################
### code chunk number 5: systemPipeR.Rnw:84-87
###################################################
library(systemPipeR)
targetspath <- paste0(system.file("extdata", package="systemPipeR"), "/targets.txt")
read.delim(targetspath, comment.char = "#")


###################################################
### code chunk number 6: systemPipeR.Rnw:91-94
###################################################
library(systemPipeR)
targetspath <- paste0(system.file("extdata", package="systemPipeR"), "/targetsPE.txt")
read.delim(targetspath, comment.char = "#")[1:2,1:6]


###################################################
### code chunk number 7: systemPipeR.Rnw:99-100
###################################################
readComp(file=targetspath, format="vector", delim="-")


###################################################
### code chunk number 8: systemPipeR.Rnw:105-107
###################################################
parampath <- paste0(system.file("extdata", package="systemPipeR"), "/tophat.param")
read.delim(parampath, comment.char = "#")


###################################################
### code chunk number 9: systemPipeR.Rnw:110-112
###################################################
args <- systemArgs(sysma=parampath, mytargets=targetspath)
args


###################################################
### code chunk number 10: systemPipeR.Rnw:115-120
###################################################
names(args)
modules(args)
cores(args)
outpaths(args)[1]
sysargs(args)[1]


###################################################
### code chunk number 11: systemPipeR.Rnw:123-124
###################################################
systemArgs(sysma=parampath, mytargets=targetspath, type="json")


###################################################
### code chunk number 12: systemPipeR.Rnw:130-131 (eval = FALSE)
###################################################
## library(systemPipeR)


###################################################
### code chunk number 13: systemPipeR.Rnw:135-137 (eval = FALSE)
###################################################
## targets <- read.delim(targetspath, comment.char = "#")
## write.table(targets, "targets_run.txt", row.names=FALSE, quote=FALSE, sep="\t")


###################################################
### code chunk number 14: systemPipeR.Rnw:141-142 (eval = FALSE)
###################################################
## args <- systemArgs(sysma=parampath, mytargets="targets_run.txt")


###################################################
### code chunk number 15: systemPipeR.Rnw:147-148 (eval = FALSE)
###################################################
## system("bowtie2-build ./data/mygenome.fa ./data/bowtie2index/mygenome")


###################################################
### code chunk number 16: systemPipeR.Rnw:152-153 (eval = FALSE)
###################################################
## bampaths <- runCommandline(args=args)


###################################################
### code chunk number 17: systemPipeR.Rnw:157-159 (eval = FALSE)
###################################################
## qsubargs <- getQsubargs(queue="batch", Nnodes="nodes=1", cores=cores(args), memory="mem=10gb", time="walltime=20:00:00")
## (joblist <- qsubRun(args=args, qsubargs=qsubargs, Nqsubs=4, package="systemPipeR"))


###################################################
### code chunk number 18: systemPipeR.Rnw:163-166 (eval = FALSE)
###################################################
## read_statsDF <- alignStats(args, fqgz=TRUE) 
## read_statsDF <- cbind(read_statsDF[targets$FileName,], targets)
## write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")


###################################################
### code chunk number 19: systemPipeR.Rnw:170-173 (eval = FALSE)
###################################################
## symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"), 
##             urlbase="http://myserver.edu/~username/", 
## 	    urlfile="IGVurl.txt")


###################################################
### code chunk number 20: systemPipeR.Rnw:178-181 (eval = FALSE)
###################################################
## parampath <- paste0(system.file("extdata", package="systemPipeR"), "/bowtieSE.param")
## args <- systemArgs(sysma=parampath, mytargets="targets_run.txt")
## bampaths <- runCommandline(args=args)


###################################################
### code chunk number 21: systemPipeR.Rnw:185-187 (eval = FALSE)
###################################################
## qsubargs <- getQsubargs(queue="batch", Nnodes="nodes=1", cores=cores(args), memory="mem=10gb", time="walltime=20:00:00")
## (joblist <- qsubRun(args=args, qsubargs=qsubargs, Nqsubs=4, package="systemPipeR"))


###################################################
### code chunk number 22: systemPipeR.Rnw:192-194 (eval = FALSE)
###################################################
## txdb <- makeTranscriptDbFromGFF(file="data/mygenome.gtf", format="gtf", dataSource="ENSEMBL", species="My_species")
## saveDb(txdb, file="./data/My_species.sqlite")


###################################################
### code chunk number 23: systemPipeR.Rnw:198-210 (eval = FALSE)
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
### code chunk number 24: systemPipeR.Rnw:215-224 (eval = FALSE)
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
### code chunk number 25: systemPipeR.Rnw:228-234 (eval = FALSE)
###################################################
## library(ape)
## rpkmDFeByg <- read.table("./results/rpkmDFeByg.xls", check.names=FALSE)
## rpkmDFeByg <- rpkmDFeByg[rowMeans(rpkmDFeByg) > 50,]
## d <- cor(rpkmDFeByg, method="spearman")
## hc <- hclust(as.dist(1-d))
## plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)


###################################################
### code chunk number 26: systemPipeR.Rnw:238-241
###################################################
targetspath <- paste0(system.file("extdata", package="systemPipeR"), "/targets.txt")
targets <- read.delim(targetspath, comment.char = "#")
(cmp <- readComp(file=targetspath, format="matrix", delim="-"))


###################################################
### code chunk number 27: systemPipeR.Rnw:243-244 (eval = FALSE)
###################################################
## edgeDF <- run_edgeR(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=TRUE, mdsplot="")


###################################################
### code chunk number 28: sessionInfo
###################################################
toLatex(sessionInfo())


