systemPipeR
===

R package to run command-line software for NG sequence analysis. Supports both
interactive job submissions and batch submissions to queuing systems of
clusters (currently tested only with Torque).

### Usage
Instructions for running _systemPipeR_ are given in the
[_vignette_](https://github.com/tgirke/systemPipeR/blob/master/vignettes/systemPipeR.pdf?raw=true).
The expected format to define the RNA-Seq samples (e.g. FASTQ files) and their
labels are given in
[_targets.txt_](https://github.com/tgirke/systemPipeR/blob/master/inst/extdata/targets.txt)
and
[_targetsPE.txt_](https://github.com/tgirke/systemPipeR/blob/master/inst/extdata/targetsPE.txt)
(latter is for PE reads). 
The run parameters of command-line software are defined by param files that have a JSON-like 
name/value structure. Here is a sample param file for _Tophat2_: [_tophat.param_](https://github.com/tgirke/systemPipeR/blob/master/inst/extdata/tophat.param). 
 
