##########################################
## Helper Functions not Exposed to User ##
##########################################

####################################################################
## Subset Reads by Mapping Regions to Create Mini FASTQ/BAM Files ##
####################################################################
## (a) E.g. download from NCBI's SRA fastq files for study SRP010938
## (b) Align reads with systemPipeR/tophat against truncated TAIR10 reference (each chr truncated to 100kbp)
## (c) Extract 100,000 reads from each fastq file mapping to the truncated regions in reference
## Step (c) is performed with the following function
.subsetReadsByMappingRegion <- function(args) {
	chromosomelength <- 100000
	mydir <- getwd()
	setwd("./data/SRP010938_sub") # outdir
	for(i in seq(along=outpaths(args))) {
		fl <- outpaths(args)[i]
		si <- seqinfo(BamFile(fl))                                                                                                                                                                                                                 
		gr <- GRanges(seqnames(si), IRanges(100, seqlengths(si)-100))                                                                                                                                                                                
		aligns <- readGAlignments(fl, param=ScanBamParam(which=gr), use.names=TRUE)
		keepids <- names(aligns[start(aligns) < chromosomelength]) # Return read ids mapping in first 100000 nucleotides of chromosomes
		myN <- sample(90000:100000, 1) # Keep number of random sampled reads between 90-100K
		keepids <- sample(unique(keepids), myN) # random sample x reads
		reads1 <- readFastq(infile1(args)[i]) # Reads in a FASTQ file from the current folder
		index <- gsub(" .*", "", as.vector(id(reads1))) %in% keepids
		reads1 <- reads1[index] # subset by keepids
		writeFastq(reads1, gsub("^.*/", "", infile1(args)[i]), full=TRUE) # writes ShortReadQ object to output file
		rm(reads1); gc() # Clean up memory
		reads2 <- readFastq(infile2(args)[i]) 
		index <- gsub(" .*", "", as.vector(id(reads2))) %in% keepids
		reads2 <- reads2[index] 
		writeFastq(reads2, gsub("^.*/", "", infile2(args)[i]), full=TRUE) 
		rm(reads2); gc() # Clean up memory
		print(i)
	}
	setwd(mydir)
}
## Usage:
# args <- systemArgs(sysma="tophat.param", mytargets="./data/SRP010938/targets.txt")
# .subsetReadsByMappingRegion(args=args)
