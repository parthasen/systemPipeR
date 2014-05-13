###############################################################
## Functions to Run NGS Aligners on Cluster or Interactively ##
###############################################################

########################################
## Specify Arguments for NGS Aligners ##
########################################
systemArgs <- function(app="tophat2", mymodules, mydir, myargs, myref, mygff, mytargets, myindir="/data/", myoutdir="/results/") {
	## Preprocessing of targets input
	mytargets <- read.delim(paste(mydir, "/", mytargets, sep=""), comment.char = "#")
	colnames(mytargets)[1] <- "FileName1" # To support FileName column for SE data
	## Insert empty FileName2 column if not present
	if(length(mytargets$FileName2)==0) mytargets <- data.frame(FileName1=mytargets$FileName1, FileName2="", mytargets[,!colnames(mytargets) %in% "FileName1"])
	
	## Tophat2
	if(app=="tophat2") {
		tophatargs <- list(modules = mymodules,
        		           args = myargs,
					# sample args: c(software="tophat", p="-p 4", g="-g 1", segment_length="--segment-length 25", i="-i 30", I="-I 3000")
                			# -G: supply GFF with transcript model info (preferred!)
					# -g: ignore all alignments with >g matches
					# -p: number of threads to use for alignment step
					# -i/-I: min/max intron lengths (50, 500000 are defaults)
					# --segment-length: length of split reads (25 is default)
                 	   	   reference = myref, 
                 	   	   gff = paste("-G ", mydir, myindir, mygff, sep=""), # assign empty string to 'mygff' if GFF/GTF is not needed
                 	   	   outpath = paste(mydir, myoutdir, sep=""), 
                 	   	   infile1 = as.character(mytargets$FileName1),
		 	   	   infile2 = as.character(mytargets$FileName2)
		 		)
		if(nchar(mygff)==0) tophatargs[["gff"]] <- "" # removes "-G" if GFF/GTF is not needed
		return(tophatargs)
	}
	
	## Bowtie2
	if(app=="bowtie2") {
		bowtie2args <- list(modules = mymodules,
        		           args = myargs,
				        # sample args: c(software="bowtie2" p="-p 4" k="-k 50", other="--non-deterministic", ")
                			# -a: report all alignments for each read
					# -k: report at most k alignments for each read 
					# --non-deterministic: more approporiat for samples with many identical reads
					# -p: number of threads to use for alignment step
				   reference = myref, 
			           outpath = paste(mydir, myoutdir, sep=""), 
                 	   	   infile1 = as.character(mytargets$FileName1),
		 	   	   infile2 = as.character(mytargets$FileName2)
		 		)
		return(bowtie2args)
	}
}
## Usage:
# mymodules <- c("bowtie2/2.1.0", "tophat/2.0.8b")
# myargs <- c(software="tophat", p="-p 4", g="-g 1", segment_length="--segment-length 25", i="-i 30", I="-I 3000")
# tophatargs <- systemArgs(app="tophat2", mymodules=mymodules, mydir=getwd(), myargs=myargs, myref="./data/My_genome.fasta", mygff="My_species.gff", mytargets="targets_run.txt", myindir="/data/", myoutdir="/results/")
# myargs <- c(software="bowtie2", p="-p 4", k="-k 50", other="--non-deterministic")
# bowtieargs <- systemArgs(app="bowtie2", mymodules=mymodules, mydir=getwd(), myargs=myargs, myref="./data/My_genome.fasta", mytargets="targets_run.txt", myindir="/data/", myoutdir="/results/")

#########################################################################
## Function to run Tophat2 including sorting and indexing of BAM files ##
#########################################################################
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
## How to run in interactive session, e.g. via qsub -I
# runTophat(tophatargs=tophatargs, runid="01")

######################################################################### 
## Function to run Bowtie2 including sorting and indexing of BAM files ##
#########################################################################
runBowtie <- function(bowtieargs=bowtieargs, runid="01") {
	library(modules); library(Rsamtools)
	moduleload(bowtieargs$modules[1]) # loads bowtie2 from module system
	tmp <- paste(bowtieargs$outpath, gsub("^.*/", "", bowtieargs$infile1), ".bam", sep="")
	bam_paths <- rep(FALSE, length(tmp)); names(bam_paths) <- tmp
	for(i in seq(along=bowtieargs$infile1)) {
		## Run alignmets only for samples for which no BAM file is available.
		bamexists <- names(bam_paths)[i]
		bam_paths[i] <- file.exists(bamexists)
		if(as.logical(bam_paths[i])) {
			next()
		} else {
			## SE and PE require different execution syntax in Bowtie2
			if(nchar(bowtieargs$infile2[i]) == 0) { 
         			bowtie_command <- paste(paste(bowtieargs$args, collapse=" "), " -x ", bowtieargs$reference, " -U ", bowtieargs$infile1[i], " -S ", bowtieargs$infile1[i], ".sam", sep="")
			} else {
         			bowtie_command <- paste(paste(bowtieargs$args, collapse=" "), " -x ", bowtieargs$reference, " -1 ", bowtieargs$infile1[i], " -2 ", bowtieargs$infile2[i], " -S ", bowtieargs$infile1[i], ".sam", sep="")
			}
			## Create submitrunID_log file
			cat(bowtie_command, file=paste(bowtieargs$outpath, "/submitargs", runid, "_log", sep=""), sep = "\n", append=TRUE)
        		## Run bowtie 
			system(bowtie_command)
			asBam(file=paste(bowtieargs$outpath, gsub("^.*/", "", bowtieargs$infile1[i]), ".sam", sep=""), destination=paste(bowtieargs$outpath, gsub("^.*/", "", bowtieargs$infile1[i]), sep=""), overwrite=TRUE, indexDestination=TRUE)
			unlink(paste(bowtieargs$outpath, gsub("^.*/", "", bowtieargs$infile1[i]), ".sam", sep=""))
		}
	}
	cat("Missing alignment results (bam files):", sum(!bam_paths), "\n"); cat("Existing alignment results (bam files):", sum(bam_paths), "\n")
	bam_paths[1:length(bam_paths)] <- file.exists(names(bam_paths))
	return(bam_paths)
}
## How to run in interactive session, e.g. via qsub -I
# runBowtie(bowtieargs=bowtieargs, runid="01")

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
# qsubargs <- getQsubargs(queue="batch", Nnodes="nodes=1", cores=as.numeric(gsub("^.* ", "", tophatargs$args["p"])), memory="mem=10gb", time="walltime=20:00:00")

##########################################################################
## Function to submit to job to the cluster (e.g. runTophat or similar) ##
##########################################################################
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
# qsubRun(appfct="runBowtie(appargs, runid)", appargs=bowtieargs, qsubargs=qsubargs, Nqsubs=1, submitdir="results", myfct="systemPipe.R")

#####################
## Alignment Stats ##
#####################
alignStats <- function(fqpaths, bampaths, fqgz=TRUE) {
	library(ShortRead); library(Rsamtools)
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
returnRPKM <- function(counts, gffsub) {
        geneLengthsInKB <- sum(width(reduce(gffsub)))/1000 # Length of exon union per gene in kbp
        millionsMapped <- sum(counts)/1e+06 # Factor for converting to million of mapped reads.
        rpm <- counts/millionsMapped # RPK: reads per kilobase of exon model.
        rpkm <- rpm/geneLengthsInKB # RPKM: reads per kilobase of exon model per million mapped reads.
        return(rpkm)
}

## Usage:
# countDFrpkm <- apply(countDF, 2, function(x) returnRPKM(counts=x, gffsub=eByg))

###############################################
## Read Sample Comparisons from Targets File ##
###############################################
## Parses sample comparisons from <COMP> line in targets.txt file. All possible
## comparisons can be specified with 'CMPset: ALL'.
readComp <- function(myfile, format="vector", delim="-") {
	if(!format %in% c("vector", "matrix")) stop("Argument format can only be vector or matrix!")
	## Parse <CMP> line
	comp <- readLines(myfile)
	comp <- comp[grepl("<CMP>", comp)]
	comp <- gsub("#.*<CMP>| {1,}", "", comp)
	comp <- strsplit(comp, ":|,")
	names(comp) <- lapply(seq(along=comp), function(x) comp[[x]][1])	
	comp <- sapply(names(comp), function(x) comp[[x]][-1], simplify=FALSE)	
	
	## Check whether all samples are present in Factor column of targets file
	checkvalues <- unique(unlist(strsplit(unlist(comp), "-")))
	checkvalues <- checkvalues[checkvalues!="ALL"]
	all <- unique(as.character(read.delim(myfile, comment.char = "#")$Factor))
	if(any(!checkvalues %in% all)) stop(paste("The following samples are not present in Factor column of targets file:", paste(checkvalues[!checkvalues %in% all], collapse=", ")))	

	## Generate outputs 
	allindex <- sapply(names(comp), function(x) any(grepl("ALL", comp[[x]])))
	if(any(allindex)) for(i in which(allindex)) comp[[i]] <- combn(all, m=2, FUN=paste, collapse=delim)
	if(format == "vector" & delim != "-") comp <- sapply(names(comp), function(x) gsub("-", delim, comp[[x]]), simplify=FALSE)
	if(format == "vector") return(comp)
	if(format == "matrix") return(sapply(names(comp), function(x) do.call("rbind", strsplit(comp[[x]], "-")), simplify=FALSE))
}
## Usage:
# cmp <- readComp(myfile="targets.txt", format="vector", delim="-")
