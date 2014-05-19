### R code from vignette source 'systemPipeR.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: systemPipeR.Rnw:41-43
###################################################
options(width=95)
unlink("test.db")


###################################################
### code chunk number 3: systemPipeR.Rnw:63-65 (eval = FALSE)
###################################################
## # $ R CMD build systemPipeR # Builds package
## install.packages("systemPipeR.1.0.0.tar.gz", repos=NULL, type="source") # Installs the package


###################################################
### code chunk number 4: systemPipeR.Rnw:70-71
###################################################
library("systemPipeR") # Loads the package


###################################################
### code chunk number 5: systemPipeR.Rnw:73-75 (eval = FALSE)
###################################################
## library(help="systemPipeR") # Lists all functions and classes 
## vignette("systemPipeR") # Opens this PDF manual from R


###################################################
### code chunk number 6: sessionInfo
###################################################
toLatex(sessionInfo())


