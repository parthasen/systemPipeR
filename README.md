systemPipeR
===

Convenience utilities to run command-line software from within R. Supports
interactive job submissions on a single machine or batch submissions to a
queuing system (currently tested only with Torque) of a compute cluster.

### Usage
Instructions for running _systemPipe.R_ are given in _analysis.R_. The expected
format to define the RNA-Seq samples (e.g. FASTQ files) and their labels
is given in _targets.txt_. Currently, supported command-line software includes:

 - Bowtie 2
 - TopHat 2 
 
