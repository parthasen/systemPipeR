### R code from vignette source 'systemPipeRNAseq.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex(use.unsrturl=FALSE)


###################################################
### code chunk number 2: systemPipeRNAseq.Rnw:43-45
###################################################
options(width=95)
unlink("test.db")


###################################################
### code chunk number 3: systemPipeRNAseq.Rnw:72-73
###################################################
library(systemPipeR)


###################################################
### code chunk number 4: systemPipeRNAseq.Rnw:77-78 (eval = FALSE)
###################################################
## source("systemPipeRNAseq_Fct.R")


###################################################
### code chunk number 5: systemPipeRNAseq.Rnw:83-86
###################################################
targetspath <- paste0(system.file("extdata", package="systemPipeR"), "/targets.txt")
targets <- read.delim(targetspath, comment.char = "#")[,1:4]
targets


###################################################
### code chunk number 6: systemPipeRNAseq.Rnw:98-103 (eval = FALSE)
###################################################
## args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
## fqlist <- seeFastq(fastq=infile1(args), batchsize=100000, klength=8)
## pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
## seeFastqPlot(fqlist)
## dev.off()


###################################################
### code chunk number 7: systemPipeRNAseq.Rnw:115-117 (eval = FALSE)
###################################################
## args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
## sysargs(args)[1] # Command-line parameters for first FASTQ file


###################################################
### code chunk number 8: systemPipeRNAseq.Rnw:120-124 (eval = FALSE)
###################################################
## moduleload(modules(args))
## system("bowtie2-build ./data/aedes-aegypti-liverpool_scaffolds_AaegL3.fa ./data/aedes-aegypti-liverpool_scaffolds_AaegL3.fa")
## qsubargs <- getQsubargs(queue="batch", cores=cores(args), memory="mem=10gb", time="walltime=20:00:00")
## (joblist <- qsubRun(args=args, qsubargs=qsubargs, Nqsubs=18, package="systemPipeR"))


###################################################
### code chunk number 9: systemPipeRNAseq.Rnw:127-128 (eval = FALSE)
###################################################
## file.exists(outpaths(args))


###################################################
### code chunk number 10: systemPipeRNAseq.Rnw:133-135 (eval = FALSE)
###################################################
## read_statsDF <- alignStats(args=args, fqgz=TRUE) 
## write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")


###################################################
### code chunk number 11: systemPipeRNAseq.Rnw:137-138 (eval = FALSE)
###################################################
## read.delim("results/alignStats.xls")


###################################################
### code chunk number 12: systemPipeRNAseq.Rnw:143-146 (eval = FALSE)
###################################################
## symLink2bam(sysargs=args, htmldir=c("~/.html/", "projects/AlexRaikhel/2014/"), 
##             urlbase="http://biocluster.ucr.edu/~tgirke/", 
## 	    urlfile="./results/IGVurl.txt")


###################################################
### code chunk number 13: systemPipeRNAseq.Rnw:152-166 (eval = FALSE)
###################################################
## library("GenomicFeatures"); library(BiocParallel)
## txdb <- loadDb("./data/AedesAegypti.sqlite")
## eByg <- exonsBy(txdb, by=c("gene"))
## bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
## multicoreParam <- MulticoreParam(workers=8); register(multicoreParam); registered()
## counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union", 
##                                                ignore.strand=TRUE, 
##                                                inter.feature=FALSE, 
##                                                singleEnd=TRUE)) 
## countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
## rownames(countDFeByg) <- names(rowData(counteByg[[1]])); colnames(countDFeByg) <- names(bfl)
## rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
## write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
## write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")


###################################################
### code chunk number 14: systemPipeRNAseq.Rnw:169-170 (eval = FALSE)
###################################################
## read.delim("results/countDFeByg.xls", row.names=1, check.names=FALSE)[1:4,1:5]


###################################################
### code chunk number 15: systemPipeRNAseq.Rnw:173-174 (eval = FALSE)
###################################################
## read.delim("results/rpkmDFeByg.xls", row.names=1, check.names=FALSE)[1:4,1:4]


###################################################
### code chunk number 16: systemPipeRNAseq.Rnw:179-187 (eval = FALSE)
###################################################
## library(ape)
## rpkmDFeByg <- read.delim("./results/rpkmDFeByg.xls", row.names=1, check.names=FALSE)[,-19]
## rpkmDFeByg <- rpkmDFeByg[rowMeans(rpkmDFeByg) > 50,]
## d <- cor(rpkmDFeByg, method="spearman")
## hc <- hclust(as.dist(1-d))
## pdf("results/sample_tree.pdf")
## plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
## dev.off()


###################################################
### code chunk number 17: systemPipeRNAseq.Rnw:198-203 (eval = FALSE)
###################################################
## library(edgeR)
## countDF <- read.delim("countDFeByg.xls", row.names=1, check.names=FALSE) 
## targets <- read.delim("targets.txt", comment="#")
## cmp <- readComp(file="targets.txt", format="matrix", delim="-")
## edgeDF <- run_edgeR(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=FALSE, mdsplot="")


###################################################
### code chunk number 18: systemPipeRNAseq.Rnw:207-212 (eval = FALSE)
###################################################
## desc <- read.delim("data/desc.xls") 
## desc <- desc[!duplicated(desc[,1]),]
## descv <- as.character(desc[,2]); names(descv) <- as.character(desc[,1])
## edgeDF <- data.frame(edgeDF, Desc=descv[rownames(edgeDF)], check.names=FALSE)
## write.table(edgeDF, "./results/edgeRglm_allcomp.xls", quote=FALSE, sep="\t", col.names = NA)


###################################################
### code chunk number 19: systemPipeRNAseq.Rnw:216-221 (eval = FALSE)
###################################################
## edgeDF <- read.delim("results/edgeRglm_allcomp.xls", row.names=1, check.names=FALSE) 
## pdf("results/DEGcounts.pdf")
## DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=1))
## dev.off()
## write.table(DEG_list$Summary, "./results/DEGcounts.xls", quote=FALSE, sep="\t", row.names=FALSE)


###################################################
### code chunk number 20: systemPipeRNAseq.Rnw:232-239 (eval = FALSE)
###################################################
## library("biomaRt") 
## listMarts() # Choose BioMart databases, here vb_mart_22 (VectorBase)
## vb <- useMart("vb_mart_22"); listDatasets(vb) # Choose genome from VectorBase, here aaegypti_eg_gene
## vb <- useMart("vb_mart_22", dataset="aaegypti_eg_gene")
## listAttributes(vb) # Choose data types you want to download
## go <- getBM(attributes=c("ensembl_gene_id",  "go_accession", "go_name_1006", "go_namespace_1003"), mart=vb)
## go[1:4,]


###################################################
### code chunk number 21: systemPipeRNAseq.Rnw:243-244 (eval = FALSE)
###################################################
## downloadGOdata(rerun=FALSE) # Do only once


###################################################
### code chunk number 22: systemPipeRNAseq.Rnw:248-262 (eval = FALSE)
###################################################
## edgeDF <- read.delim("results/edgeRglm_allcomp.xls", row.names=1, check.names=FALSE) 
## DEGlist <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=5))
## up_down <- DEGlist$UporDown; names(up_down) <- paste(names(up_down), "_up_down", sep="")  
## up <- DEGlist$Up; names(up) <- paste(names(up), "_up", sep="")  
## down <- DEGlist$Down; names(down) <- paste(names(down), "_down", sep="")  
## DEGlist <- c(up_down, up, down)
## DEGlist <- DEGlist[sapply(DEGlist, length) > 0]
## loadData("data/GO")
## BatchResult <- GOCluster_Report(setlist=DEGlist, method="all", id_type="gene", CLSZ=10, cutoff=0.9, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
## write.table(BatchResult, "./results/GOBatchResultedgeR_allcomp.xls", quote=FALSE, sep="\t", col.names = NA)
## library("biomaRt"); vb <- useMart("vb_mart_22", dataset="aaegypti_eg_gene")
## goslimvec <- as.character(getBM(attributes=c("goslim_goa_accession"), mart=vb)[,1])
## BatchResultslim <- GOCluster_Report(setlist=DEGlist, method="slim", id_type="gene", myslimv=goslimvec, CLSZ=10, cutoff=0.01, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
## write.table(BatchResultslim, "./results/GOslimBatchResultedgeR_allcomp.xls", quote=FALSE, sep="\t", col.names = NA)


###################################################
### code chunk number 23: systemPipeRNAseq.Rnw:266-271 (eval = FALSE)
###################################################
## gos <- read.delim("results/GOslimBatchResultedgeR_allcomp.xls", row.names=1, check.names=FALSE)
## gos <- gos[grep("^iEcRa24h-iLuc24h", gos$CLID), ] 
## pdf("./results/GOslimbarplotMF.pdf", height=8, width=10); goBarplot(gos, gocat="MF"); dev.off()
## pdf("./results/GOslimbarplotBP.pdf", height=8, width=10); goBarplot(gos, gocat="BP"); dev.off()
## pdf("./results/GOslimbarplotCC.pdf", height=8, width=10); goBarplot(gos, gocat="CC"); dev.off()


###################################################
### code chunk number 24: sessionInfo
###################################################
toLatex(sessionInfo())


