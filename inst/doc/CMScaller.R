## ----prepareSession, include=FALSE--------------------------------------------
library(Biobase)
library(BiocStyle)
knitr::opts_chunk$set(fig.width=6, fig.height=3, 
        dev.args=list(pointsize=8), dpi=150,
        collapse=TRUE, message=TRUE, echo=TRUE, warnings=FALSE)
options(scipen=-1, digits=2)

## ---- eval=FALSE--------------------------------------------------------------
#  # dependencies: run if not already installed
#  source("https://bioconductor.org/biocLite.R")
#  biocLite(c("Biobase", "limma"))
#  # proper repository to be fixed for publication
#  install.packages("pathToPackageFile/CMScaller_0.99.2.tar.gz", repos=NULL)

## ----quickStart, fig.cap="CMScaller graphic output. Left heatmap shows relative expression for template genes. Samples (columns) are ordered according to class predictions and confidence. The height of the white bars below gives the unadjusted prediction $p$-values. Genes (rows) are ordered according to class template. Heatmap color saturation indicates magnitude while red and blue indicate genes up and down relative to the sample mean. Right heatmap shows results for Camera gene set analysis. Heatmap color saturation indicates statistical significance and red and blue indicates direction of change."----
library(Biobase) # if input is ExpressionSet
library(CMScaller)
# get RNA-seq counts from TCGA example data
counts <- exprs(crcTCGAsubset)
head(counts[,1:2])
# prediction and gene set analysis
par(mfrow=c(1,2))
res <- CMScaller(emat=counts, RNAseq=TRUE, FDR=0.05)
cam <- CMSgsa(emat=counts, class=res$prediction,RNAseq=TRUE)
# comparison with true class
table(pred=res$prediction, true=crcTCGAsubset$CMS)
head(res, n=3)

## ----makeTemplates, fig.keep="last", fig.height=4-----------------------------
emat <- crcTCGAsubset
cms <- emat$CMS.Syn
train <- sample(seq_along(cms), size=length(cms)/(2))
deg <- subDEG(emat[,train], class=cms[train], doVoom=TRUE)
templates <- ntpMakeTemplates(deg, resDEG=TRUE, topN=50)
templates$symbol <- fromTo(templates$probe)
tail(templates,n=3)

## ----visGSA, message=TRUE, fig.cap="Gene Set Analysis (GSA) shows that CMS are biologically distinct.", fig.width=3----
# increase left margins to accommodate gene set names
par.old <- par()
par(mfrow=c(1,1), mar=par.old$mar+c(0,4,0,0))
subCamera(counts, cms, geneList=geneSets.CRC, doVoom=TRUE)
# restore margins
par(mar=par.old$mar)

## ----visPCA, message=TRUE, fig.cap="Principal component analysis (PCA) and CMS. First two principal components seperates CMS (left) with CMS4 characgterized by high levels of THBS4 and low levels of CLCA1 (right)."----
# increase left margins to accommodate gene set names
par(mfrow=c(1,2))
p <- subPCA(emat = crcTCGAsubset, class = crcTCGAsubset$CMS.Syn, 
            normMethod = "quantile", pch=16, frame=FALSE)
plotPC(p, n=6, entrez=TRUE)

## ----input--------------------------------------------------------------------
# loads included emat, scales and centers
emat <- crcTCGAsubset
emat_sc <- ematAdjust(emat, normMethod="quantile")
head(emat_sc[,1:2])

## -----------------------------------------------------------------------------
# test set prediction
res <- ntp(emat_sc[,-train], templates, nPerm=1000)
res <- subSetNA(res, pValue=.1)
table(pred=res$prediction, true=cms[-train])
head(res)

## ----principleDistance, fig.width=3-------------------------------------------
# random centered/scaled expression matrix and templates
set.seed(42)
N <- 5000;P <- 5000;K <- 4;nPerm <- 1000;n <- 1
X <- matrix(rnorm(P*N, mean=0, sd=1), ncol=N)
Y <- matrix(rbinom(P*K, size=1, prob=.01), ncol=K)
# sample-template correlations (implemented in corCosine)
cos.sim <- crossprod(X,Y) / outer(
                sqrt(apply(X, 2, crossprod)), 
                sqrt(apply(Y, 2, crossprod)))
# sample-template distances (vectorized)
simToDist <- function(cos.sim) sqrt(1/2 * (1-cos.sim))
cos.dist <- simToDist(cos.sim)
hist(cos.dist, xlab="cosine correlation distance")

## ----principlePermutations----------------------------------------------------
# estimate prediction confidence
pred.class <- apply(cos.dist, 1, which.min)
pred.dists <- apply(cos.dist, 1, min)
null.dist <- replicate(nPerm, min(simToDist(corCosine(sample(X[,n]), Y))))
p <- rank(c(pred.dists[n], null.dist))[1]/(length(null.dist))

## ----pUniform, fig.cap="NTP results for random data. Left heatmap shows expression for template genes for random data with rows and columns sorted according to class. Right histogram shows the expected uniform $p$-value distribution for the random data."----
# rearrange matrix and templates for ntp input
rownames(X) <- make.names(seq_len(P))
templates <- lapply(seq_len(K), function(k) rownames(X)[Y[,k]==1])
names(templates) <- paste("k", seq_len(K))
templates <- ntpMakeTemplates(templates, resDEG=FALSE)
# permutations set to 100 to reduce processing for vignette
par(mfrow=c(1,2))
res <- ntp(X, templates, nCores=1L, nPerm=100, doPlot=TRUE)
# expect uniform distribution
hist(res$p.value, main ="", xlab="prediction p-values")

## ----endSession, echo=FALSE, results='asis'-----------------------------------
toLatex(sessionInfo())

