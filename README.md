**`CMScaller`: an R package for consensus molecular subtyping of colorectal cancer pre-clinical models**

The `CMScaller` package provides Consensus Molecular Subtype (CMS) classification of colorectal cancer pre-clinical models [Guinney 2015; Eide 2017]. A small ensembl of functions for evaluating and visualizing results is also included. The core algorithm is *[Nearest Template Prediction (NTP)](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0015543)* algorithm proposed by Yujin Hoshida [Hoshida 2010]. See accompagnying vignette for further details.

**Install**

From R/RStudio uncomment as required and run the following code 

```{r}
### dependencies: run if not already installed
### limma has lof of dependencies - takes time
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("Biobase", "limma"))
# install.packages("devtools", repos = NULL)

### install: latest version
devtools::install_github("Lothelab/CMScaller")
```

**Quick start**

```{r}
library(Biobase)
library(CMScaller)
par(mfrow=c(1,2))

### CMS prediction of TCGA primary colorectal cancers
res <- CMScaller(crcTCGAsubset, RNAseq=TRUE, doPlot=TRUE)
head(res)

### Camera Gene Set Analysis with CMS informative gene sets
cam <- CMSgsa(emat=crcTCGAsubset, class=res$prediction, RNAseq=TRUE)
head(cam$CMS4)

### limma differential gene expression analysis and visualization
deg <- subDEG(emat=crcTCGAsubset, class=res$prediction, doVoom=TRUE)
subVolcano(deg, geneID="symbol")
```

**Design**

Package builds on *[Bioconductor](http://bioconductor.org/)* to facilitate integration with existing workflows [Huber 2015] and was developed in *[RStudio](https://www.rstudio.com/)* following guidelines in [R packages](http://r-pkgs.had.co.nz/) [Wickham 2015]. Devtools and roxygen

**References**

Eide PW, Bruun J, Lothe RA, Sveen A. CMScaller: stroma-independent consensus molecular subtyping of colorectal cancer pre-clinical models. Sci Rep. 2017.

Guinney J, Dienstmann R, Wang X, de Reynies A, Schlicker A, Soneson C, et al. The consensus molecular subtypes of colorectal cancer. Nat Med. 2015;21:1350-6.

Hoshida, Y. Nearest Template Prediction: A Single-Sample-Based Flexible Class Prediction with Confidence Assessment. PLoS ONE 2010;5, e15543.

Huber W, Carey VJ, Gentleman R, Anders S, Carlson M, Carvalho BS, et al. Orchestrating high-throughput genomic analysis with Bioconductor. Nat Meth. 2015;12:115–21. 

R Core Team. R: A Language and Environment for Statistical Computing [Internet]. Vienna, Austria: R Foundation for Statistical Computing; 2013. Available from: http://www.R-project.org/

Wickham, H. R Packages: Organize, Test, Document, and Share Your Code. 1st ed. O’Reilly Media. 2015.

Wickham, H. Chang, W. devtools: Tools to Make Developing R Packages Easier. R package version 1.13.0. 2017.

Wickham, H., Danenberg, P., Eugster, M. roxygen2: In-Line Documentation for R. R package version 6.0.1. 2017.

