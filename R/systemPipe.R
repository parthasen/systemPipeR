###############################################################
## Functions to Run NGS Aligners on Cluster or Interactively ##
###############################################################

##############################################
## Class and Method Definitions for SYSargs ##
##############################################
## Define SYSargs class
setClass("SYSargs", representation(modules="character", 
				   software="character",
				   cores="numeric", 
				   other="character",
				   reference="character",
			           infile1="character",
				   infile2="character",
				   outfile1="character",
				   sysargs="character", 
				   outpaths="character")
)

## Methods to return SYSargs as list
setGeneric(name="modules", def=function(x) standardGeneric("modules"))
setMethod(f="modules", signature="SYSargs", definition=function(x) {return(as.character(x@modules))})
setGeneric(name="software", def=function(x) standardGeneric("software"))
setMethod(f="software", signature="SYSargs", definition=function(x) {return(as.character(x@software))})
setGeneric(name="cores", def=function(x) standardGeneric("cores"))
setMethod(f="cores", signature="SYSargs", definition=function(x) {return(x@cores)})
setGeneric(name="other", def=function(x) standardGeneric("other"))
setMethod(f="other", signature="SYSargs", definition=function(x) {return(x@other)})
setGeneric(name="reference", def=function(x) standardGeneric("reference"))
setMethod(f="reference", signature="SYSargs", definition=function(x) {return(x@reference)})
setGeneric(name="infile1", def=function(x) standardGeneric("infile1"))
setMethod(f="infile1", signature="SYSargs", definition=function(x) {return(x@infile1)})
setGeneric(name="infile2", def=function(x) standardGeneric("infile2"))
setMethod(f="infile2", signature="SYSargs", definition=function(x) {return(x@infile2)})
setGeneric(name="outfile1", def=function(x) standardGeneric("outfile1"))
setMethod(f="outfile1", signature="SYSargs", definition=function(x) {return(x@outfile1)})
setGeneric(name="SampleName", def=function(x) standardGeneric("SampleName"))
setMethod(f="SampleName", signature="SYSargs", definition=function(x) {return(names(x@sysargs))})
setGeneric(name="sysargs", def=function(x) standardGeneric("sysargs"))
setMethod(f="sysargs", signature="SYSargs", definition=function(x) {return(x@sysargs)})
setGeneric(name="outpaths", def=function(x) standardGeneric("outpaths"))
setMethod(f="outpaths", signature="SYSargs", definition=function(x) {return(x@outpaths)})

## Constructor methods
## List to SYSargs with: as(mylist, "SYSargs")
setAs(from="list", to="SYSargs",  
        def=function(from) {
		new("SYSargs", modules=from$modules, 
			       software=from$software,
			       cores=from$cores, 
			       other=from$other,
			       reference=from$reference,
			       infile1=from$infile1, 
			       infile2=from$infile2,
			       outfile1=from$outfile1,
			       sysargs=from$sysargs, 
                               outpaths=from$outpaths)
})

## Define print behavior for SYSargs
setMethod(f="show", signature="SYSargs", 
	definition=function(object) {    
	cat("An instance of ", class(object), " with ", length(object@sysargs), " samples ", "\n", sep="")
})

## Extend names() method
setMethod(f="names", signature="SYSargs",
    definition=function(x) {
    	return(slotNames(x))
})

## Construct SYSargs object from param and targets files
systemArgs <- function(sysma, mytargets, type="SYSargs") {
	## Read sysma and convert to arglist
	sysma <- as.matrix(read.delim(sysma, comment.char = "#"))
	sysma[is.na(sysma)] <- ""
	arglist <- sapply(as.character(unique(sysma[,"PairSet"])), function(x) as.vector(t(as.matrix(sysma[sysma[,"PairSet"]==x, 2:3]))))
	for(i in seq(along=arglist)) names(arglist[[i]]) <- paste(rep(c("n", "v"), length(arglist[[i]])/2), rep(1:(length(arglist[[i]])/2), 2), sep="")
	if(type=="json") return(toJSON(arglist))
	## Validity checks
	mytargets <- read.delim(mytargets, comment.char = "#")
	if(any(duplicated(mytargets$SampleName))) stop("SampleName column of mytargets cannot contain duplicated entries!")
	## Preprocessing of targets input
	colnames(mytargets)[1] <- "FileName1" # To support FileName column for SE data
	## Insert empty FileName2 column if not present
	if(length(mytargets$FileName2)==0) mytargets <- data.frame(FileName1=mytargets$FileName1, FileName2="", mytargets[,!colnames(mytargets) %in% "FileName1"])
	## Check name:value violations in arglist
	check <- sapply(names(arglist), function(x) sum(grepl("^n", names(arglist[[x]]))) == sum(grepl("^n", names(arglist[[x]]))))
	if(any(!check)) stop(paste("Name:Value violation in arglist component(s):", paste(names(check[check]), collapse=", ")))
	
	## Modify arglist object as specified in arglist and mytargets
	## Remove module component and store values in separate container
	modules <- as.character(arglist$modules[grepl("v", names(arglist$modules))])
	arglist <- arglist[!names(arglist) %in% "modules"]
	
	## Extract single value components
	software <- as.character(arglist$software[grepl("v", names(arglist$software))])
	other <- as.character(arglist$other[grepl("v", names(arglist$other))])
	reference <- as.character(arglist$reference[grepl("v", names(arglist$reference))])
	cores <- as.numeric(arglist$cores[grepl("v", names(arglist$cores))])	

	## Populate arglist$infile1
	infile1 <- gsub("<|>", "", arglist$infile1[grepl("^<.*>$", arglist$infile1)][[1]])
	infile1 <- as.character(mytargets[,infile1])	
	argname <- arglist$infile1[grep("<.*>", arglist$infile1)[1] -1]
	path <- arglist$infile1[grep("path", arglist$infile1)[1] +1]
	infile1back <- paste(path, infile1, sep="")
	names(infile1back) <- as.character(mytargets$SampleName)	
	infile1 <- paste(argname, " ", path, infile1, sep="")
	arglist[["infile1"]] <- gsub("(^ {1,})|( ${1,})", "", infile1)
	
	## Populate arglist$infile2
	infile2 <- gsub("<|>", "", arglist$infile2[grepl("^<.*>$", arglist$infile2)][[1]])
	infile2 <- as.character(mytargets[,infile2])	
	argname <- arglist$infile2[grep("<.*>", arglist$infile2)[1] -1]
	path <- arglist$infile2[grep("path", arglist$infile2)[1] +1]
	infile2back <- paste(path, infile2, sep="")
	names(infile2back) <- as.character(mytargets$SampleName)	
	infile2 <- paste(argname, " ", path, infile2, sep="")
	arglist[["infile2"]] <- gsub("(^ {1,})|( ${1,})", "", infile2)
	
	## Populate arglist$outfile1
	outfile1 <- gsub("<|>", "", arglist$outfile1[grepl("^<.*>$", arglist$outfile1)][[1]])
	outfile1 <- as.character(mytargets[,outfile1])
	outfile1 <- gsub("^.*/", "", outfile1) 	
	remove <- arglist$outfile1[grep("remove", arglist$outfile1)[1] +1]
	outfile1 <- gsub(as.character(remove), "", outfile1)
	outfile1back <- outfile1
	outpaths <- outfile1
	outextension <- as.character(arglist$outfile1[grep("outextension", arglist$outfile1)+1])
	append <- arglist$outfile1[grep("append", arglist$outfile1)[1] +1]
	outfile1 <- paste(outfile1, append, sep="")
	argname <- arglist$outfile1[grep("<.*>", arglist$outfile1)[1] -1]
	path <- arglist$outfile1[grep("path", arglist$outfile)[1] +1]
	outfile1back <- paste(path, outfile1, sep="")
	names(outfile1back) <- as.character(mytargets$SampleName)	
	outfile1 <- paste(argname, " ", path, outfile1, sep="")
	arglist[["outfile1"]] <- gsub("(^ {1,})|( ${1,})", "", outfile1)

	## Generate arglist$outpaths
	outpaths <- paste(path, outpaths, outextension, sep="")
	names(outpaths) <- as.character(mytargets$SampleName)	

	## Collapse remaining components to single string vectors
	remaining <- names(arglist)[!names(arglist) %in% c("outfile1", "infile1", "infile2", "outpaths")]
	for(i in remaining) arglist[[i]] <- rep(gsub("(^ {1,})|( ${1,})", "", paste(arglist[[i]], collapse=" ")), length(arglist$infile1))
	args <- do.call("cbind", arglist)	
	rownames(args) <- as.character(mytargets$SampleName)
	args <- apply(args, 1, paste, collapse=" ")
	
	## Construct SYSargs object from components
	syslist <- list(modules=modules, 
                        software=software, 
                        cores=cores,
			other=other,
			reference=reference,
			infile1=infile1back,
			infile2=infile2back,
			outfile1=outfile1back,
                        sysargs=args, 
                        outpaths=outpaths)
	sysargs <- as(syslist, "SYSargs")
	if(type=="SYSargs") return(sysargs)
}
## Usage:
# args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
# names(args); modules(args); cores(args); outpaths(args); sysargs(args)

##############################################################################
## Function to run NGS aligners including sorting and indexing of BAM files ##
##############################################################################
runCommandline <- function(args, runid="01") {
	for(j in modules(args)) moduleload(j) # loads specified software from module system
	commands <- sysargs(args)
	completed <- file.exists(outpaths(args))
	names(completed) <- outfile1(args)
	resultdir <- gsub("(\\./.*?/).*", "\\1", outpaths(args)[[1]])	
	for(i in seq(along=commands)) {
		## Run alignmets only for samples for which no BAM file is available.
		if(as.logical(completed)[i]) {
			next()
		} else {
			## Create submitrunID_log file
			cat(commands[i], file=paste(resultdir, "submitargs", runid, "_log", sep=""), sep = "\n", append=TRUE)
        		## Run executable  
			system(as.character(commands[i]))
			if(grepl(".bam$", names(completed[i]))) { # If output is unindexed *.bam file (e.g. Tophat2)
				sortBam(file=names(completed[i]), destination=gsub("\\.bam$", "", names(completed[i])))
        			indexBam(names(completed[i]))
			}
			if(grepl(".sam$", names(completed[i]))) { # If output is *.sam file (e.g. Bowtie2)
				asBam(file=names(completed[i]), destination=gsub("\\.sam$", "", names(completed[i])), overwrite=TRUE, indexDestination=TRUE)
				unlink(names(completed[i]))
			}
		}
	}
	bamcompleted <- gsub("sam$", "bam$", file.exists(outpaths(args)))
	names(bamcompleted) <- SampleName(args)
	cat("Missing alignment results (bam files):", sum(!as.logical(bamcompleted)), "\n"); cat("Existing alignment results (bam files):", sum(as.logical(bamcompleted)), "\n")
	return(bamcompleted)
}

## Usage: 
# runCommandline(args=args)

####################
## qsub Arguments ##
####################
getQsubargs <- function(software="qsub", queue="batch", Nnodes="nodes=1", cores=as.numeric(gsub("^.* ", "", tophatargs$args["p"])), memory="mem=10gb", time="walltime=20:00:00") {
	qsubargs <- list(software=software, 
			queue=queue, 
                	Nnodes=Nnodes, 
                 	cores=cores, 
                 	memory=memory, 
                 	time=time)
	return(qsubargs)
}
## Usage:
# qsubargs <- getQsubargs(queue="batch", Nnodes="nodes=1", cores=cores(tophat), memory="mem=10gb", time="walltime=20:00:00")

###########################################################################
## Function to submit runCommandline jobs to queuing system of a cluster ##
###########################################################################
qsubRun <- function(appfct, appargs, qsubargs, Nqsubs=1, submitdir="results", package="systemPipeR") {
	appargs <- sysargs(appargs)
	mydir <- getwd()
	setwd(submitdir)
	splitvector <- sort(rep_len(1:Nqsubs, length.out=length(appargs)))
	commands <- split(appargs, splitvector)
	qsub_command <- paste(qsubargs$software, " -q ", qsubargs$queue, " -l ", qsubargs$Nnodes, ":ppn=", qsubargs$cores, ",", qsubargs$memory, ",", qsubargs$time, sep="")	
	jobids <- NULL
	for(i in 1:Nqsubs) {
		appargs[["infile1"]] <- filesets[[i]]
		appargs[["infile2"]] <- filesets2[[i]]
		counter <- formatC(i, width = 2, format = "d", flag = "0")
		appfct <- gsub("runid.*)", paste("runid=", "'", counter, "'", ")", sep=""), appfct) # Passes on proper runid
		save(list=c("appargs"), file=paste("submitargs", counter, sep="")) 
		rscript <- c(paste("library('", package, "')", sep=""), paste("load('submitargs", counter, "')", sep=""), appfct)
		writeLines(rscript, paste("submitargs", counter, ".R", sep=""))
		writeLines(c("#!/bin/bash", "cd $PBS_O_WORKDIR", paste("Rscript --verbose submitargs", counter, ".R", sep="")), paste("submitargs", counter, ".sh", sep=""))
		myqsub <- paste(qsub_command, paste("submitargs", counter, ".sh", sep=""))
		(jobids <- c(jobids, system(myqsub, intern=TRUE)))
	}
	setwd(mydir)
	names(filesets) <- jobids
	return(filesets)
}
## Usage:
# qsubRun(appfct="runCommandline(args=args)", appargs=tophat, qsubargs=qsubargs, Nqsubs=1, submitdir="results", package="systemPipeR")

#####################
## Alignment Stats ##
#####################
alignStats <- function(fqpaths, bampaths, fqgz=TRUE) {
	fqpaths <- fqpaths[as.logical(bampaths)]
	bampaths <- bampaths[as.logical(bampaths)]
	## Obtain total read number from FASTQ files
	if(fqgz==TRUE) {
		Nreads <- sapply(fqpaths, function(x) as.numeric(system(paste("zcat", x, "| wc -l | cut -d' ' -f1"), intern=TRUE))/4)
	} else {
		Nreads <- sapply(fqpaths, function(x) as.numeric(system(paste("wc -l", x, "| cut -d' ' -f1"), intern=TRUE))/4)
	}
	## Obtain total number of alignments from BAM files
	bfl <- BamFileList(names(bampaths), yieldSize=50000, index=character())
	Nalign <- countBam(bfl)
	## Obtain number of primary alignments from BAM files
	param <- ScanBamParam(flag=scanBamFlag(isNotPrimaryRead=FALSE, isUnmappedQuery=FALSE))
	Nalignprim <- countBam(bfl, param=param)
	statsDF <- data.frame(FileName=names(Nreads), 
                              Nreads=Nreads, 
                              Nalign=Nalign$records, 
                              Perc_Aligned=Nalign$records/Nreads*100, 
                              Nalign_Primary=Nalignprim$records, 
                              Perc_Aligned_Primary=Nalignprim$records/Nreads*100
	)
	return(statsDF)
}
## Usage:
#read_statsDF <- alignStats(fqpaths=tophatargs$infile1, bampaths=bampaths, fqgz=TRUE) 

########################
## RPKM Normalization ##
########################
returnRPKM <- function(counts, ranges) {
        geneLengthsInKB <- sum(width(reduce(ranges)))/1000 # Length of exon union per gene in kbp
        millionsMapped <- sum(counts)/1e+06 # Factor for converting to million of mapped reads.
        rpm <- counts/millionsMapped # RPK: reads per kilobase of exon model.
        rpkm <- rpm/geneLengthsInKB # RPKM: reads per kilobase of exon model per million mapped reads.
        return(rpkm)
}

## Usage:
# countDFrpkm <- apply(countDF, 2, function(x) returnRPKM(counts=x, ranges=eByg))

###############################################
## Read Sample Comparisons from Targets File ##
###############################################
## Parses sample comparisons from <CMP> line(s) in targets.txt file. All possible
## comparisons can be specified with 'CMPset: ALL'.
readComp <- function(file, format="vector", delim="-") {
	if(!format %in% c("vector", "matrix")) stop("Argument format can only be assigned: vector or matrix!")
	## Parse <CMP> line
	comp <- readLines(file)
	comp <- comp[grepl("<CMP>", comp)]
	comp <- gsub("#.*<CMP>| {1,}", "", comp)
	comp <- strsplit(comp, ":|,")
	names(comp) <- lapply(seq(along=comp), function(x) comp[[x]][1])	
	comp <- sapply(names(comp), function(x) comp[[x]][-1], simplify=FALSE)	
	
	## Check whether all samples are present in Factor column of targets file
	checkvalues <- unique(unlist(strsplit(unlist(comp), "-")))
	checkvalues <- checkvalues[checkvalues!="ALL"]
	all <- unique(as.character(read.delim(file, comment.char = "#")$Factor))
	if(any(!checkvalues %in% all)) stop(paste("The following samples are not present in Factor column of targets file:", paste(checkvalues[!checkvalues %in% all], collapse=", ")))	

	## Generate outputs 
	allindex <- sapply(names(comp), function(x) any(grepl("ALL", comp[[x]])))
	if(any(allindex)) for(i in which(allindex)) comp[[i]] <- combn(all, m=2, FUN=paste, collapse=delim)
	if(format == "vector" & delim != "-") comp <- sapply(names(comp), function(x) gsub("-", delim, comp[[x]]), simplify=FALSE)
	if(format == "vector") return(comp)
	if(format == "matrix") return(sapply(names(comp), function(x) do.call("rbind", strsplit(comp[[x]], "-")), simplify=FALSE))
}
## Usage:
# cmp <- readComp(file="targets.txt", format="vector", delim="-")

#################################
## Access module system from R ##
#################################
## List software available in module system
modulelist <- function() {
	system("bash -c \"module avail\"")
}

## Load software from module system
moduleload <- function(module) { 
	modpath <- system(paste("bash -c \"module load", module, "; export | grep '^declare -x PATH='\""), intern = TRUE) 
	modpath <- gsub("^declare -x PATH=\"(.*)\"$", "\\1", modpath, perl = TRUE)
	Sys.setenv(PATH=modpath) 
}

#######################################################################
## Run edgeR GLM with entire count matrix or subsetted by comparison ##
#######################################################################
## If independent=TRUE then countDF will be subsetted for each comparison
run_edgeR <- function(countDF, targets, cmp, independent=TRUE, paired=NULL, mdsplot="") {
    samples <- as.character(targets$Factor); names(samples) <- paste(as.character(targets$SampleName), "", sep="")
    countDF <- countDF[, names(samples)]
    countDF[is.na(countDF)] <- 0
    edgeDF <- data.frame(row.names=rownames(countDF))
    group <- as.character(samples)
    if(independent==TRUE) {
        loopv <- seq(along=cmp[,1])
    } else {
	loopv <- 1
    }
    for(j in loopv) {
	## Filtering and normalization
	y <- DGEList(counts=countDF, group=group) # Constructs DGEList object
	if(independent == TRUE) {
	    subset <- samples[samples %in% cmp[j,]]
	    y <- y[, names(subset)]
        }
	keep <- rowSums(cpm(y)>1) >= 2; y <- y[keep, ]
	y <- calcNormFactors(y)
	## Design matrix
	if(length(paired)==0) {
		design <- model.matrix(~0+y$samples$group, data=y$samples)
		colnames(design) <- levels(y$samples$group)
	} else {
        	if(length(paired)>0 & independent==FALSE) stop("When providing values under 'paired' also set independent=TRUE")
		Subject <- factor(samplepairs[samples %in% cmp[j,]])
		Treat <- y$samples$group
        	design <- model.matrix(~Subject+Treat)
		levels(design) <- levels(y$samples$group)
	}
        ## Estimate dispersion
	y <- estimateGLMCommonDisp(y, design, verbose=TRUE) # Estimates common dispersions
	y <- estimateGLMTrendedDisp(y, design) # Estimates trended dispersions
	y <- estimateGLMTagwiseDisp(y, design) # Estimates tagwise dispersions 
	fit <- glmFit(y, design) # Fits the negative binomial GLM for each tag and produces an object of class DGEGLM with some new components.
	## Contrast matrix is optional but makes anlysis more transparent
	if(independent == TRUE) {
		mycomp <- paste(cmp[j,1], cmp[j,2], sep="-")
	} else {
		mycomp <- paste(cmp[,1], cmp[,2], sep="-")
	}
	if(length(paired)==0) contrasts <- makeContrasts(contrasts=mycomp, levels=design)
	for(i in seq(along=mycomp)) {
	    if(length(paired)==0) {
            	lrt <- glmLRT(fit, contrast=contrasts[,i]) # Takes DGEGLM object and carries out the likelihood ratio test. 
	    } else {
                lrt <- glmLRT(fit) # No contrast matrix with paired design
            }
            deg <- as.data.frame(topTags(lrt, n=length(rownames(y))))
	    colnames(deg) <- paste(paste(mycomp[i], collapse="_"), colnames(deg), sep="_")
	    edgeDF <- cbind(edgeDF, deg[rownames(edgeDF),]) 
	}
	if(nchar(mdsplot)>0) pdf(paste("./results/sample_MDS_", paste(unique(subset), collapse="-"), ".pdf", sep="")); plotMDS(y); dev.off()
    }
    return(edgeDF)
}
## Usage:
# cmp <- readComp(file=targetspath, format="matrix", delim="-")
# edgeDF <- run_edgeR(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=TRUE, mdsplot="")


