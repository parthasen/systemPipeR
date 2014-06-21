systemPipeR
===

R package for building *end-to-end* analysis pipelines with automated report
generation for NGS applications such as RNA-Seq, ChIP-Seq, VAR-Seq and many
others. An important feature is support for running command-line software, such
as NGS aligners, on both single machines or compute clusters. This includes
both interactive job submissions or batch submissions to queuing systems of
clusters (tested only with Torque).

### Usage
Instructions for running _systemPipeR_ are given in the
[_vignette_](https://github.com/tgirke/systemPipeR/blob/master/vignettes/systemPipeR.pdf?raw=true).
The sample data set used in the vignette can be downloaded [_here_](http://biocluster.ucr.edu/~tgirke/projects/systemPipeR_test_data.zip). 
The expected format to define the RNA-Seq samples (e.g. FASTQ files) and their
labels are given in
[_targets.txt_](https://github.com/tgirke/systemPipeR/blob/master/inst/extdata/targets.txt)
and
[_targetsPE.txt_](https://github.com/tgirke/systemPipeR/blob/master/inst/extdata/targetsPE.txt)
(latter is for PE reads). 
The run parameters of command-line software are defined by param files that have a simplified
JSON-like name/value structure. Here is a sample param file for _Tophat2_: [_tophat.param_](https://github.com/tgirke/systemPipeR/blob/master/inst/extdata/tophat.param). 
