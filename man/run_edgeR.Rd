\name{run_edgeR}
\alias{run_edgeR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Runs edgeR 
}
\description{
Convenience function to rund \code{edgeR GLM} on many comparisons. Analysis 
can be performed on entire count matrix or only on those columns required for each
comparison. 
}
\usage{
run_edgeR(countDF, targets, cmp, independent = TRUE, paired = NULL, mdsplot = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{countDF}{
\code{date.frame} containing raw read counts
}
  \item{targets}{
	targets \code{data.frame}
}
  \item{cmp}{
	\code{character matrix} where comparisons are defined in two columns. This matrix should be generated with \code{readComp()} from targets file. Values used for comparisons need to match those in \code{Factor} column of targets file.
}
  \item{independent}{
	If \code{independent=TRUE} then \code{countDF} will be subsetted for each comparison
}
  \item{paired}{
	Defines pairs (\code{character vector}) for paired analysis. Default is unpaired (\code{paired=NULL}).
}
  \item{mdsplot}{
	Directory where \code{plotMDS} should be written to. Default setting \code{mdsplot=""} will omit the plotting step.
}
}
\details{
...
}
\value{
\code{data.frame} containing \code{edgeR} results from all comparisons. Comparison labels are appended to column titles for tracking.
}
\references{
...
}
\author{
Thomas Girke
}
\seealso{
\code{edgeR} vignette and \code{readComp}
}
\examples{
\dontrun{
cmp <- readComp(file=targetspath, format="matrix", delim="-")
edgeDF <- run_edgeR(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=TRUE, mdsplot="")
pval <- edgeDF[, grep("_FDR$", colnames(edgeDF)), drop=FALSE]
fold <- edgeDF[, grep("_logFC$", colnames(edgeDF)), drop=FALSE]
DEG_list <- filterDEGs(pval=pval, log2FC=fold, filter=c(Fold=2, FDR=1), plot=TRUE)
names(DEG_list)
DEG_list$Summary
}
}
\keyword{ utilities }
