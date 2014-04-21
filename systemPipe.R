#########################################################
## Functions to Run Tophat on cluster or interactively ##
#########################################################
## Bowtie2/Tophat2 arguments
systemArgs <- function(mymodules, mydir, myargs, myref, mygff, mytargets, myindir="/data/", myoutdir="/results/") {
	tophatargs <- list(modules = mymodules,
        	           args = myargs,
				# sample args: c(software="tophat", p="-p 4", g="-g 1", segment_length="--segment-length 25", i="-i 30", I="-I 3000")
                		# -G: supply GFF with transcript model info (preferred!)
				# -g: ignore all alginments with >g matches
				# -p: number of threads to use for alignment step
				# -i/-I: min/max intron lengths (50, 500000 are defaults)
				# --segment-length: length of split reads (25 is default)
                 	   reference = paste(mydir, myindir, myref, sep=""), 
                 	   gff = paste("-G ", mydir, myindir, mygff, sep=""), # assign empty string if GFF/GTF is not needed
                 	   outpath = paste(mydir, myoutdir, sep=""), 
                 	   infile1 = as.character(read.delim(paste(mydir, "/", mytargets, sep=""), comment.char = "#")$FileName),
		 	   infile2 = rep("", length(read.delim(paste(mydir, "/", mytargets, sep=""), comment.char = "#")$FileName))
		 	)
	if(nchar(mygff)==0) tophatargs[["gff"]] <- "" # removes "-G" if GFF/GTF is not needed
	return(tophatargs)
}
## Usage:
# mymodules <- c("bowtie2/2.1.0", "tophat/2.0.8b")
# myargs <- c(software="tophat", p="-p 4", g="-g 1", segment_length="--segment-length 25", i="-i 30", I="-I 3000")
# tophatargs <- systemArgs(mymodules=mymodules, mydir=getwd(), myargs=myargs, myref="TAIR10_chr_all.fas", mygff="TAIR10_GFF3_genes.gff", mytargets="targets_run.txt")

## Function to run Bowtie2/Tophat2 including sorting and indexing of BAM files
runTophat <- function(tophatargs=tophatargs, runid="01") {
	library(modules); library(Rsamtools)
	moduleload(tophatargs$modules[1]); moduleload(tophatargs$modules[2]) # loads bowtie2/tophat2 from module system
	tmp <- paste(tophatargs$outpath, gsub("^.*/", "", tophatargs$infile1), ".tophat/accepted_hits.bam", sep="")
	bam_paths <- rep(FALSE, length(tmp)); names(bam_paths) <- tmp
	for(i in seq(along=tophatargs$infile1)) {
		## Run alignmets only for samples for which no BAM file is available.
		bamexists <- names(bam_paths)[i]
		bam_paths[i] <- file.exists(bamexists)
		if(as.logical(bam_paths[i])) {
			next()
		} else {
        		tophat_command <- paste(paste(tophatargs$args, collapse=" "), " ", tophatargs$gff, " -o ", tophatargs$outpath, gsub("^.*/", "", tophatargs$infile1[i]), ".tophat ", tophatargs$reference, " ", tophatargs$infile1[i], " ", tophatargs$infile2[i], sep="")
			## Create submitrunID_log file
			cat(tophat_command, file=paste(tophatargs$outpath, "/submitargs", runid, "_log", sep=""), sep = "\n", append=TRUE)
        		## Run tophat 
			system(tophat_command)
			sortBam(file=paste(tophatargs$outpath, gsub("^.*/", "", tophatargs$infile1[i]), ".tophat/accepted_hits.bam", sep=""), destination=paste(tophatargs$outpath, gsub("^.*/", "", tophatargs$infile1[i]), ".tophat/accepted_hits", sep=""))
        		indexBam(paste(tophatargs$outpath, gsub("^.*/", "", tophatargs$infile1[i]), ".tophat/accepted_hits.bam", sep=""))
		}
	}
	cat("Missing alignment results (bam files):", sum(!bam_paths), "\n"); cat("Existing alignment results (bam files):", sum(bam_paths), "\n")
	bam_paths[1:length(bam_paths)] <- file.exists(names(bam_paths))
	return(bam_paths)
}
## How to run in interactive session, e.g. on owl or with qsub -I
# runTophat(tophatargs=tophatargs, runid="01")

## qsub arguments
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
# qsubargs <- getQsubargs(queue="batch", Nnodes="nodes=1", cores=as.numeric(gsub("^.* ", "", tophatargs$args["p"])), memory="mem=10gb", time="walltime=20:00:00")

## Function to submit runTophat (or similar) to cluster
qsubRun <- function(appfct="runTophat(appargs, runid)", appargs=tophatargs, qsubargs=qsubargs, Nqsubs=1, submitdir="results", myfct="systemPipe.R") {
	mydir <- getwd()
	setwd(submitdir)
	splitvector <- sort(rep_len(1:Nqsubs, length.out=length(appargs$infile1)))
	filesets <- split(appargs$infile1, splitvector)
	qsub_command <- paste(qsubargs$software, " -q ", qsubargs$queue, " -l ", qsubargs$Nnodes, ":ppn=", qsubargs$cores, ",", qsubargs$memory, ",", qsubargs$time, sep="")	
	jobids <- NULL
	for(i in 1:Nqsubs) {
		appargs[["infile1"]] <- filesets[[i]]
		counter <- formatC(i, width = 2, format = "d", flag = "0")
		appfct <- gsub("runid.*)", paste("runid=", "'", counter, "'", ")", sep=""), appfct) # Passes on proper runid
		save(list=c("appargs"), file=paste("submitargs", counter, sep="")) 
		rscript <- c(paste("source('", myfct, "')", sep=""), paste("load('submitargs", counter, "')", sep=""), appfct)
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
# qsubRun(appfct="runTophat(appargs, runid)", appargs=tophatargs, qsubargs=qsubargs, Nqsubs=1, submitdir="results", myfct="systemPipe.R")

## Alignment Stats
alignStats <- function(fqpaths, bampaths, fqgz=TRUE) {
	library(ShortRead); library(Rsamtools)
	fqpaths <- fqpaths[as.logical(bampaths)]
	bampaths <- bampaths[as.logical(bampaths)]
	if(fqgz==TRUE) {
		Nreads <- sapply(fqpaths, function(x) as.numeric(system(paste("zcat", x, "| wc -l | cut -d' ' -f1"), intern=TRUE))/4)
	} else {
		Nreads <- sapply(fqpaths, function(x) as.numeric(system(paste("wc -l", x, "| cut -d' ' -f1"), intern=TRUE))/4)
	}
	bfl <- BamFileList(names(bampaths), yieldSize=50000, index=character())
	Nalign <- countBam(bfl)
	statsDF <- data.frame(FileName=names(Nreads), Nreads=Nreads, Nalign=Nalign$records, Perc_Aligned=Nalign$records/Nreads*100)
	return(statsDF)
}
## Usage:
#read_statsDF <- alignStats(fqpaths=tophatargs$infile1, bampaths=bampaths, fqgz=TRUE) 

## RPKM Normalization
returnRPKM <- function(counts, gffsub) {
        geneLengthsInKB <- sum(width(reduce(gffsub)))/1000 # Length of exon union per gene in kbp
        millionsMapped <- sum(counts)/1e+06 # Factor for converting to million of mapped reads.
        rpm <- counts/millionsMapped # RPK: reads per kilobase of exon model.
        rpkm <- rpm/geneLengthsInKB # RPKM: reads per kilobase of exon model per million mapped reads.
        return(rpkm)
}

## Usage:
#countDFrpkm <- apply(countDF, 2, function(x) returnRPKM(counts=x, gffsub=eByg))

## Usage:
#countDFrpkm <- apply(countDF, 2, function(x) returnRPKM(counts=x, gffsub=eByg))

