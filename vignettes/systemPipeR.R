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
### code chunk number 13: systemPipeR.Rnw:135-136 (eval = FALSE)
###################################################
## args <- systemArgs(sysma="tophat.param", mytargets="targetsPE.txt")


###################################################
### code chunk number 14: systemPipeR.Rnw:145-149 (eval = FALSE)
###################################################
## fqlist <- seeFastq(fastq=infile1(args), batchsize=10000, klength=8)
## pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
## seeFastqPlot(fqlist)
## dev.off()


###################################################
### code chunk number 15: systemPipeR.Rnw:160-162 (eval = FALSE)
###################################################
## moduleload(modules(args)) # Skip if module system is not available
## system("bowtie2-build ./data/tair10.fasta ./data/tair10.fasta")


###################################################
### code chunk number 16: systemPipeR.Rnw:166-167 (eval = FALSE)
###################################################
## bampaths <- runCommandline(args=args)


###################################################
### code chunk number 17: systemPipeR.Rnw:171-173 (eval = FALSE)
###################################################
## qsubargs <- getQsubargs(queue="batch", cores=cores(args), memory="mem=10gb", time="walltime=20:00:00")
## (joblist <- qsubRun(args=args, qsubargs=qsubargs, Nqsubs=18, package="systemPipeR"))


###################################################
### code chunk number 18: systemPipeR.Rnw:177-179 (eval = FALSE)
###################################################
## read_statsDF <- alignStats(args, fqgz=TRUE) 
## write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")


###################################################
### code chunk number 19: systemPipeR.Rnw:183-186 (eval = FALSE)
###################################################
## symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"), 
##             urlbase="http://myserver.edu/~username/", 
## 	    urlfile="IGVurl.txt")


###################################################
### code chunk number 20: systemPipeR.Rnw:191-193 (eval = FALSE)
###################################################
## args <- systemArgs(sysma="bowtieSE.param", mytargets="targets.txt")
## bampaths <- runCommandline(args=args)


###################################################
### code chunk number 21: systemPipeR.Rnw:197-199 (eval = FALSE)
###################################################
## qsubargs <- getQsubargs(queue="batch", cores=cores(args), memory="mem=10gb", time="walltime=20:00:00")
## (joblist <- qsubRun(args=args, qsubargs=qsubargs, Nqsubs=18, package="systemPipeR"))


###################################################
### code chunk number 22: systemPipeR.Rnw:204-207 (eval = FALSE)
###################################################
## library(GenomicFeatures)
## txdb <- makeTranscriptDbFromGFF(file="data/tair10.gff", format="gff", dataSource="TAIR", species="A. thaliana")
## saveDb(txdb, file="./data/tair10.sqlite")


###################################################
### code chunk number 23: systemPipeR.Rnw:211-222 (eval = FALSE)
###################################################
## library(BiocParallel)
## txdb <- loadDb("./data/tair10.sqlite")
## eByg <- exonsBy(txdb, by="gene")
## bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
## multicoreParam <- MulticoreParam(workers=4); register(multicoreParam); registered()
## counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union", ignore.strand=TRUE, inter.feature=TRUE, singleEnd=TRUE)) # Note: for strand-specific RNA-Seq set 'ignore.strand=FALSE' and for PE data set 'singleEnd=FALSE'
## countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
## rownames(countDFeByg) <- names(rowData(counteByg[[1]])); colnames(countDFeByg) <- names(bfl)
## rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
## write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
## write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")


###################################################
### code chunk number 24: systemPipeR.Rnw:227-236 (eval = FALSE)
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
### code chunk number 25: systemPipeR.Rnw:240-246 (eval = FALSE)
###################################################
## library(ape)
## rpkmDFeByg <- read.table("./results/rpkmDFeByg.xls", check.names=FALSE)
## rpkmDFeByg <- rpkmDFeByg[rowMeans(rpkmDFeByg) > 50,]
## d <- cor(rpkmDFeByg, method="spearman")
## hc <- hclust(as.dist(1-d))
## plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)


###################################################
### code chunk number 26: systemPipeR.Rnw:256-260
###################################################
library(edgeR)
targets <- read.delim(targetspath, comment="#")
cmp <- readComp(file=targetspath, format="matrix", delim="-")
cmp[[1]]


###################################################
### code chunk number 27: systemPipeR.Rnw:263-264 (eval = FALSE)
###################################################
## edgeDF <- run_edgeR(countDF=countDFeByg, targets=targets, cmp=cmp[[1]], independent=FALSE, mdsplot="")


###################################################
### code chunk number 28: systemPipeR.Rnw:267-268 (eval = FALSE)
###################################################
## DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=10))


###################################################
### code chunk number 29: systemPipeR.Rnw:276-278 (eval = FALSE)
###################################################
## names(DEG_list)
## DEG_list$Summary


###################################################
### code chunk number 30: sessionInfo
###################################################
toLatex(sessionInfo())


