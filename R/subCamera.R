#' limma/voom camera gene set analysis visualization
#' @export
#' @description Wrapper function to perform \code{\link[limma]{camera}} Gene
#' Set Analysis (GSA) and visualize results in barplot (two classes) or
#' heatmap (>2 classes).
#' @param emat numeric matrix with row features and sample columns.
#' @param class a factor vector specifying sample classes
#' \code{length(class)==ncol(emat)}.
#' @param batch a factor vector specifying additional level in design matrix.
#' @param keepN a logical or numeric vector specifying which samples to keep
#' (defaults to all).
#' @param geneList a named list. where each item is a gene set consisting of a
#' character vector of genes belonging  to that set.
#' @param xKey a character vector matching emat with geneList identifiers.
#' \code{length(xKey)==nrow(emat)}. Leave blank if rownames(emat) and geneList
#' uses same identifiers.
#' @param pValue a p-value (number) indicating minimum gene set significance
#' included in plot.
#' @param topN integer, number of gene sets to include in the plot.
#' (for K>2, lowest p-value for each class is included).
#' @param allowNegCor logical, passed to limma::\code{\link[limma]{camera}} as
#' allow.neg.cor parameter.
#' @param interGeneCor a number, passed to limma::\code{\link[limma]{camera}}
#' as inter.gene.cor parameter.
#' @param doPlot a logical, indicating whether to return plot (TRUE) or
#' significance values (FALSE) of results.
#' @param doVoom a logical, indicating whether emat is \emph{untransformed}
#' sequencing count data and \code{\link[limma]{voom}} should be used for
#' modeling.
#' @param normMethod a character, only used if doVoom=TRUE and passed to
#' \code{\link[edgeR]{calcNormFactors}} if element in c("TMM","RLE",
#' "upperquartile","none") or \code{\link[limma]{voom}} ("scale", "quantile",
#' "cyclicloess").
#' @param pMax a numeric, \eqn{log10(-pValue)} for color scale.
#' @param rowCluster logical, indicating whether heatmap rows should be
#' clustered.
#' @param classCol a character vector specifying class colors.
#' @param heatCol a character vector specifying heatmap colors.
#' @param cexText numeric, text scaling factor.
#' @param legendAdd a logical, whether to add legend to barplot.
#' @param ... additional arguments passed to \code{\link[graphics]{image}} or
#' in case of \code{levels(class)==2} \code{\link[graphics]{barplot}}.
#' @return a barplot/heatmap and \code{\link[limma]{camera}} output (list,
#' invisible). For two-classes comparison, gene sets are ranked by significance
#' and bars are colored according to relative up-regulation. In heatmap with
#' default \code{heatCol}, red and blue indicates relative up- and
#' down-regulation respectively. Color intensity reflects significance.
#' Nominal `camera` \eqn{p}-values are used as input for visualization.
#' @references Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucl. Acids Res. 2015;gkv007.
#' @references Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics. 2010;26:139-40.
#' @references Law CW, Chen Y, Shi W, Smyth GK. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology. 2014;15:R29.
#' @references Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. Nucleic Acids Res. 2012;gks461.
#' @examples
#' # sample subset for reduced run-time
#' subset <- 25:75
#' cam <- subCamera(emat=crcTCGAsubset, class=crcTCGAsubset$CMS,
#'     rowCluster=FALSE, doVoom=TRUE, keepN=subset)
subCamera <- function(emat, class, keepN = TRUE, batch = NULL,
                        geneList = NULL, xKey = NULL,
                        pValue = .01, topN=15,
                        allowNegCor = FALSE, interGeneCor = .01,
                        doPlot = TRUE, doVoom = FALSE,
                        normMethod = "quantile",
                        pMax = 10, rowCluster = TRUE,
                        classCol = getOption("subClassCol"),
                        heatCol = NULL, cexText = 1,
                        legendAdd = TRUE, ...)
    {

    if (!packageExists("limma"))
        stop (paste0("Function requires limma package available from:\n",
                    "http://www.bioconductor.org/packages/release/bioc/html/limma.html"))

    textfun <- function(..., cex.text)
        graphics::text(..., xpd=TRUE, cex=cexText*cex.text)

    requireNamespace("limma", quietly = TRUE)
    methods.edgeR <- c("TMM", "RLE", "upperquartile", "none")
    methods.voom <- c("scale", "quantile", "cyclicloess")

    # checkInput ##############################################################

    fallback.xKey <- NULL
    if (class(emat)[1] == "ExpressionSet") {
        fallback.xKey <- suppressPackageStartupMessages(Biobase::fData(emat))
        emat <- Biobase::exprs(emat)
    }

    cc <- stats::complete.cases(emat)
    emat <- emat[cc,]

    # imports Reactome gene sets if no geneList provided
    if (is.null(geneList)) {
        utils::data("geneSets.CMS", envir = environment())
        geneList <- CMScaller::geneSets.CMS
    }

    if(is.null(xKey)) {
        index <- lapply(geneList, function(x) rownames(emat) %in% x)
        # if suspicous matching, try to guess correct fData(emat) column
        if (sum(index[[1]]) < 3 & !is.null(fallback.xKey)) {
            guess <- apply(fallback.xKey, 2, function(x)
                sum(x %in% as.character(geneList[[1]])))
            message (paste(colnames(fallback.xKey)[which.max(guess)]),
                    " guessed input id for matching emat and geneList")
            xKey <- as.character(fallback.xKey[,which.max(guess)])
            index <- lapply(geneList, function(x) xKey %in% x)
        }
    } else {
        index <- lapply(geneList, function(x) as.character(xKey)[cc] %in% x)
    }

    if (sum(sapply(index, sum)) == 0)
        stop("unable to match emat and geneList, check identifiers")

    # drop samples with unknown class label
    if (sum(is.na(class[keepN])) > 0 ) {
        message(paste0(sum(is.na(class[keepN])),
                       " samples with class or batch NA's excluded"))
        if (is.numeric(keepN)) keepN <- seq_len(ncol(emat)) %in% keepN
        if (!is.null(batch)) {
            keepN <- keepN & !is.na(class) & !is.na(batch)
        } else {
            keepN <- keepN & !is.na(class)
        }
    }

    # coerce syntactically valid class labels
    non.empty.level <- table(class[keepN], useNA="ifany")>0
    class <- class[keepN]
    class <- factor(make.names(class))
    if (!is.null(batch)) batch <- factor(make.names(batch))

    # modelDEG ################################################################

    subtypes <- levels(class)
    emat <- emat[,keepN]
    K <- length(subtypes)
    cont.cam <- vector("list", length = length(subtypes))
    names(cont.cam) <- subtypes

    # nullBatch ###############################################################
    if (is.null(batch)) {
        design <- stats::model.matrix(~0+class, data = class)
        colnames(design) <- levels(class)

        # modelVoom ###########################################################

        if (doVoom == TRUE) emat <- voomTransform(emat, normMethod, design)

        # performCamera #######################################################
        # do group g against all other K-k
        for(k in seq_len(K)) {
            compare <- paste0(subtypes[k], "-(",
                            paste(subtypes[-k], collapse = "+"), ")/",
                            K-1)
            contrasts2 <- limma::makeContrasts(contrasts = compare, levels=design)
            cont.cam[[k]]  <- limma::camera(emat,index = index,
                                    design = design, contrast = contrasts2,
                                    inter.gene.cor = interGeneCor,
                                    sort = FALSE)
        }

    # !nullBatch ##############################################################
    } else {
        batch <- batch[keepN]
        df <- droplevels(data.frame(group = class, batch = batch))
        design <- stats::model.matrix(~0+group+batch, data = df)
        colnames(design) <- gsub("group", "", colnames(design))

        if (doVoom == TRUE) emat <- voomTransform(emat, normMethod, design)
        mm <- levels(df$group)

        # performCamera #######################################################
        # do group g against all other G-g
        for(k in seq_len(K)) {
            contrasts1 <- paste0(mm[k],"-(", paste0(mm[-k],
                                    collapse = "+"),")/",length(mm)-1)
            contrasts2 <- limma::makeContrasts(contrasts =
                                        contrasts1, levels = colnames(design))
            cont.cam[[k]]  <- limma::camera(emat,index = index,
                                        design = design, contrast = contrasts2,
                                        inter.gene.cor = interGeneCor,
                                        sort = FALSE)
        }
    }

    # prepareMatrix ###########################################################
    cam.mat <- matrix(nrow=length(index), ncol = K)

    if (K == 2) {
        cam.mat <- cont.cam[[1]] # symmetric doesn't matter which one we choose
        keep <- which(cam.mat$PValue < pValue)
        keepCam <- nrow(cam.mat[keep,])
    } else {
        for (k in seq_len(K)){
            cam.mat[,k] <- -log10(cont.cam[[k]]$PValue)
            isDown <- cont.cam[[k]]$Direction == "Down"
            cam.mat[isDown, k] <- -(-log10(cont.cam[[k]]$PValue[isDown]))
    }

    cam.mat[is.na(cam.mat)] <- 0
    colnames(cam.mat) <- subtypes
    rownames(cam.mat) <- names(index)

    keepCam <- apply(cam.mat, 1, function(x)
        max(abs(x),na.rm=TRUE)) > -log10(pValue)
    # get the per-class top hits
    if (!is.null(topN) & topN < nrow(cam.mat)) {
        cam.rank <- apply(-cam.mat, 2, rank)
        idx <- sapply(seq_len(topN), function(i) {
            unique(as.vector(apply(cam.rank, 2,function(x) which(x< i))))
            })
        idx <- idx[[which(sapply(idx, length) >= topN)[1]]]
        keepCam <- seq_len(nrow(cam.mat)) %in% idx
    } else {
        keepCam <- rep(TRUE, nrow(cam.mat))
    }

    # prepareHeatmap ##########################################################
    if( sum(keepCam) < 2) {
        message("no significant gene sets identified")
        doPlot <- FALSE
    }

    }

    # widen right margin for long gene set names
    if (isTRUE(doPlot)) {
        user.par <- graphics::par()
        on.exit(graphics::par(mar=user.par$mar))
        graphics::par(mar=graphics::par()$mar+c(0,5,0,0))
    }

    if (doPlot == TRUE & K > 2) {

        # get colors
        if (is.null(classCol)) classCol <- subData[["classCol"]]
        if (is.null(heatCol)) heatCol <- subData[["hmCol"]]

        # color breaks
        breaksLeft <- seq(-pMax,-.1, length = length(heatCol)/2)
        breaks <- c(breaksLeft, 0, abs(rev(breaksLeft)))
        labCol <- paste0(colnames(cam.mat), "\nN=", table(class))

        # plotting subset
        m <- cam.mat[keepCam,]
        P <- nrow(m)
        if (isTRUE(rowCluster)) m <- m[stats::hclust(stats::dist((m)))$order,]

        # correct break overflow (white if not fixed)
        m[m > pMax] <- pMax
        m[m < -pMax] <- -pMax

        # image cell sizes
        yy <- seq(0, 1, length.out=nrow(m)+1)
        xx <- seq(0, 1, length.out = K+1)

        # heatmap
        graphics::image(y=yy, x = xx, z=t(m),
                        col = heatCol, breaks=breaks,
                        xlim=c(0,1), ylim=c(0,1),
                        xaxt="n", yaxt="n",
                        xlab="", ylab="", ...)

        # add row and column labels
        yy <- seq(0+(1/P/2), 1-(1/P/2), length.out=nrow(m))
        textfun(-.05, yy, rownames(m), pos=2, adj=1, offset=0, cex.text=1)

        xx <- seq(0+(1/K/2), 1-(1/K/2), length.out = K)
        textfun(xx, -.05, subtypes, pos=2, offset=0.25,adj=c(.5,.5), srt=90, cex.text=1)
        xxmid <- xx

        # add classColor bar
        yy <- line2user(0:1, side=3)
        yy <- seq(yy[1], yy[2],length=3)[-1]
        xx <- seq(0, 1, length.out = K+1)
        graphics::rect(xleft=xx[-(K+1)], xright = xx[-1],
            ybottom = yy[1], ytop = yy[2], border = FALSE,
            col=classCol[non.empty.level], xpd=TRUE, lwd=.75)

        # add samples/group
        textfun(xxmid, yy[1], cex.text=.75, pos=3, offset=0,
                c(as.vector(table(class))))

        # add colorScale legend
        xx <- seq(line2user(6, 2),line2user(2, 2), length.out=3)
        yy <- rev(line2user(2:1,1))
        yy <- seq(yy[1], yy[2], length=3)
        bb <- length(breaks)

        graphics::rect(
            xleft = seq(xx[1], xx[3], length.out=bb)[-bb],
            xright = seq(xx[1], xx[3], length.out=bb)[-1],
            ybottom = yy[2], ytop=yy[3], col=heatCol, border=NA, xpd=TRUE)

        textfun(xx, (yy[2]+yy[3])/2, c(pMax,0,pMax), pos=1, cex.text=.75)
        textfun(xx[-2], yy[3], c("dn", "up"), pos=c(2,4), cex.text=.75)
        textfun(xx[2], line2user(3,1), expression(-log[10](italic(p))),
                pos=1, cex.text=.9)
    }

    if (doPlot == TRUE & K == 2) {
        if (is.null(classCol)) classCol <- subData[["classCol"]]
        if (is.null(heatCol)) heatCol <- subData[["hmCol"]]
        heatCol <- heatCol[c(1,length(heatCol))]

        if (!is.null(topN)) keep <- which(rank(cam.mat$PValue) <= topN)

        cam <- cam.mat[keep,]
        cam <- cam[order(cam$PValue, decreasing = TRUE, na.last = TRUE),]
        cam$col <- classCol[ifelse(cam$Direction == "Up", 1, 2)]

        bp <- graphics::barplot(-log10(cam$PValue), horiz = TRUE,
                names.arg = rownames(cam), cex.names = cexText,
                las = 1,
                xlab = expression(-log[10](p)),
                col = cam$col, ...)
        if (isTRUE(legendAdd)) {
            graphics::legend(0,max(bp)+1, xjust=.5,yjust=.5,
                            legend = levels(class), xpd=TRUE,horiz=TRUE,
                            fill = classCol, bty = "n", cex=cexText*.75)
        }
    }

    # returnMatrix ############################################################
    invisible(cont.cam)
}
