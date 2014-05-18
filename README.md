systemPipeR
===

R package to run command-line software for NG sequence analysis 
on both single machines or compute clusters. Supports interactive 
job submissions or batch submissions to queuing systems of clusters
(currently tested only with Torque).

### Usage
Instructions for running [_systemPipeR_](https://github.com/tgirke/systemPipeR/blob/master/vignettes/analysis.R)
are given in _./vignettes/analysis.R_. The expected
format to define the RNA-Seq samples (e.g. FASTQ files) and their labels
are given in _./inst/extdata/targets.txt_ and _./inst/extdata/targetsPE.txt_ (latter is for PE reads). Currently, 
supported command-line software includes:

 - Bowtie 2
 - TopHat 2 
 
