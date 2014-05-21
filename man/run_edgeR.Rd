\name{run_edgeR}
\alias{run_edgeR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Runs edgeR 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
run_edgeR(countDF = countDF, targets = targets, cmp = cmp, independent = TRUE, paired = NULL, mdsplot = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{countDF}{
%%     ~~Describe \code{countDF} here~~
}
  \item{targets}{
%%     ~~Describe \code{targets} here~~
}
  \item{cmp}{
%%     ~~Describe \code{cmp} here~~
}
  \item{independent}{
%%     ~~Describe \code{independent} here~~
}
  \item{paired}{
%%     ~~Describe \code{paired} here~~
}
  \item{mdsplot}{
%%     ~~Describe \code{mdsplot} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (countDF = countDF, targets = targets, cmp = cmp, independent = TRUE, 
    paired = NULL, mdsplot = "") 
{
    samples <- as.character(targets$Factor)
    names(samples) <- paste(as.character(targets$SampleName), 
        "", sep = "")
    countDF <- countDF[, names(samples)]
    countDF[is.na(countDF)] <- 0
    edgeDF <- data.frame(row.names = rownames(countDF))
    group <- as.character(samples)
    if (independent == TRUE) {
        loopv <- seq(along = cmp[, 1])
    }
    else {
        loopv <- 1
    }
    for (j in loopv) {
        y <- DGEList(counts = countDF, group = group)
        if (independent == TRUE) {
            subset <- samples[samples \%in\% cmp[j, ]]
            y <- y[, names(subset)]
        }
        keep <- rowSums(cpm(y) > 1) >= 2
        y <- y[keep, ]
        y <- calcNormFactors(y)
        if (length(paired) == 0) {
            design <- model.matrix(~0 + y$samples$group, data = y$samples)
            colnames(design) <- levels(y$samples$group)
        }
        else {
            if (length(paired) > 0 & independent == FALSE) 
                stop("When providing values under 'paired' also set independent=TRUE")
            Subject <- factor(samplepairs[samples \%in\% cmp[j, 
                ]])
            Treat <- y$samples$group
            design <- model.matrix(~Subject + Treat)
            levels(design) <- levels(y$samples$group)
        }
        y <- estimateGLMCommonDisp(y, design, verbose = TRUE)
        y <- estimateGLMTrendedDisp(y, design)
        y <- estimateGLMTagwiseDisp(y, design)
        fit <- glmFit(y, design)
        if (independent == TRUE) {
            mycomp <- paste(cmp[j, 1], cmp[j, 2], sep = "-")
        }
        else {
            mycomp <- paste(cmp[, 1], cmp[, 2], sep = "-")
        }
        if (length(paired) == 0) 
            contrasts <- makeContrasts(contrasts = mycomp, levels = design)
        for (i in seq(along = mycomp)) {
            if (length(paired) == 0) {
                lrt <- glmLRT(fit, contrast = contrasts[, i])
            }
            else {
                lrt <- glmLRT(fit)
            }
            deg <- as.data.frame(topTags(lrt, n = length(rownames(y))))
            colnames(deg) <- paste(paste(mycomp[i], collapse = "_"), 
                colnames(deg), sep = "_")
            edgeDF <- cbind(edgeDF, deg[rownames(edgeDF), ])
        }
        if (nchar(mdsplot) > 0) 
            pdf(paste("./results/sample_MDS_", paste(unique(subset), 
                collapse = "-"), ".pdf", sep = ""))
        plotMDS(y)
        dev.off()
    }
    return(edgeDF)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
