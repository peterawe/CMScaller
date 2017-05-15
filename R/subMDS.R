#' MDS plot with class labels
#' @export
#' @description Multi-dimensional scaling (MDS) plot with class labels. Wrapper
#' for \code{\link[limma]{plotMDS}}.
#' @param emat a numeric matrix with row features and sample columns.
#' @param class a factor vector specifying subtypes compared
#' \code{length(class)==ncol(emat)}.
#' @param keepN a logical or numeric vector specifying which samples to keep
#' (defaults to all).
#' @param normMethod a character, passed to \code{\link[edgeR]{calcNormFactors}}
#' if element in c("TMM","RLE", "upperquartile","none") or
#' \code{\link[limma]{voom}} ("scale", "quantile", "cyclicloess").
#' @param nGenes an integer specifying number of features to used for MDS
#' @param dim a numeric vector of length 2, specifying which principal
#' components to plot.
#' @param labelSamp a logical indicating whether points (FALSE) or
#' colnames(emat) are plotted.
#' @param labelCenters a logical indicating whether class centers should be
#' labeled.
#' @param classConfusion a factor vector with same length and levels as
#' \code{class} specifying alternative classifications.
#' @param classCol a character vector of hexcolors with
#' \code{length(classCol)>=levels(class)}. If NULL, a pre-defined RColorBrewer
#' based palette is used.
#' @param ... arguments to be passed to \code{\link[limma]{plotMDS}}.
#' @details \code{plotMDS} provides a MDS plot of emat with points colored
#' according to class labels. If \emph{classConfusion} is specified, samples
#' where \code{class != classConfusion} are highlighted.
#' @return A labeled MDS plot
#' @examples
#' # sample subset for reduced run-time
#' subset <- 25:50
#' subMDS(emat=crcTCGAsubset, class=crcTCGAsubset$CMS,
#'         keepN=subset, normMethod="quantile")
subMDS <- function(emat, class = NULL, keepN = TRUE,
                    nGenes = 1000,
                    normMethod = NULL,
                    dim = c(1,2),
                    labelSamp = FALSE,
                    labelCenters = TRUE, classConfusion = NULL,
                    classCol = getOption("subClassCol"),...)
    {

    if (labelSamp == TRUE) labelCenters <- FALSE
    if (class(emat) == "ExpressionSet") emat <- Biobase::exprs(emat)

    # prepare input data
    if (!is.null(class)) {
        class <- droplevels(as.factor(class)[keepN])
        pchCol <- classCol[class]
        pchCol[is.na(pchCol)] <- "#000000"
    } else {pchCol="#000000"}

    emat <- emat[,keepN]

    # MDS
    if (!is.null(normMethod)) emat <- voomTransform(emat, normMethod)
    if (isTRUE(labelSamp)) {
        p <- limma::plotMDS(emat, col=pchCol, dim.plot=dim,
                        top=nGenes,...)
    } else {
        colnames(emat) <- NULL
        p <- limma::plotMDS(emat, col=pchCol, dim.plot=dim,
                        top=nGenes,...)
    }

    if (!is.null(class)) {
        graphics::legend("topleft", legend = levels(class), bg = "gray95",
           col = classCol, cex = .75, pch = "*", text.col = classCol)
    }


    # add error labels
    if (!is.null(classConfusion) & !is.null(class)) {
        pscores <- cbind(p$x, p$y)
        classConfusion <- classConfusion[keepN]
        highl <- class != classConfusion
        highl[is.na(highl) | is.na(class)] <- FALSE

        dotdotdot <- list(...)
        if (!is.null(dotdotdot$cex)) cex.fac <- dotdotdot$cex else cex.fac=1

        graphics::points(pscores[highl,1],pscores[highl,2],
                         pch=23, bg=0, cex=cex.fac)
        graphics::points(pscores[highl,1],pscores[highl,2],
                         pch=23, cex=cex.fac,
                         bg=paste0(classCol[class[highl]], "3F"))
        graphics::points(pscores[highl,1],pscores[highl,2],
                         pch=23, cex=cex.fac,
                         bg=paste0(classCol[classConfusion[highl]], "3F"))
        graphics::legend("bottomright",  pch=23, pt.bg="gray", cex =.75,
                         legend="discordance")
    }
}
