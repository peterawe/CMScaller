#' PCA plot with class labels.
#' @description Principal component analysis (PCA) plot with class labels.
#' @export
#' @param emat numeric matrix with row features and sample columns.
#' @param class a factor vector specifying sample classes
#' \code{length(class)==ncol(emat)}.
#' @param keepN a logical or numeric vector specifying which samples to keep.
#' @param nGenes an integer specifying number of features to used for PCA.
#' @param dim a numeric vector of length 2, specifying which principal
#' components to plot.
#' @param normMethod a character, passed to \code{\link[edgeR]{calcNormFactors}}
#' if element in c("TMM","RLE", "upperquartile","none") or
#' \code{\link[limma]{voom}} ("scale", "quantile", "cyclicloess").#'
#' @param labelSamp a logical indicating whether points or
#' \code{colnames(emat)} are plotted.
#' @param legend a character specifying legend placement (e.g. "topleft",
#' "bottomright").
#' @param classConfusion a factor vector with same length and levels as
#' \code{class} specifying alternative classifications.
#' @param classCol a character vector of hexcolors with
#' \code{length(classCol)>=levels(class)}.
#' @param ... arguments to be passed to \code{\link[graphics]{plot}}.
#' @details \code{plotPCA} provides a PCA with points colored
#' according to class labels. If \code{classConfusion} is specified, samples
#' where \code{class != classConfusion} are highlighted.
#' @return a labeled PCA plot and list of class \code{prcomp} passed
#' from \code{\link[stats]{prcomp}}
#' @seealso \code{\link[stats]{prcomp}}, \code{\link{subMDS}}
subPCA <- function(emat, class = NULL, keepN = TRUE,
                    nGenes = 1000,
                    normMethod = NULL,
                    dim = c(1,2),
                    labelSamp = FALSE, legend = "topright",
                    classConfusion = NULL,
                    classCol = getOption("subClassCol"),...)
    {

    if (class(emat) == "ExpressionSet") emat <- Biobase::exprs(emat)

    # prepare input data
    if (!is.null(class)) {
        class <- droplevels(as.factor(class)[keepN])
        K <- length(levels(droplevels(class)))
        pchCol <- classCol[class]
        pchCol[is.na(pchCol)] <- "#000000"
    } else {pchCol="#000000"}

    emat <- emat[,keepN]

    if (!is.null(normMethod)) emat <- voomTransform(emat, normMethod)$E

    # PCA, IPR is the 10-90% inter-percentile-range
    pemat <-stats::prcomp((t(emat[rank(-apply(emat, 1, IPR)) <= nGenes,])))
    pscores <- pemat$x[,dim]

    # prepare for plotting
    percVar <- round(as.vector((pemat$sdev^2)/sum(pemat$sdev^2)*100),0)

    dotdotdot <- list(...)
    if (!is.null(dotdotdot$cex)) cex.fac <- dotdotdot$cex else cex.fac=1

    # plot
    if (labelSamp==FALSE) {
    graphics::plot(pscores, col = pchCol, ...,
        xlab = paste0("PC", dim[1], " ", percVar[dim[1]],"%"),
        ylab = paste0("PC", dim[2], " ", percVar[dim[2]],"%"))
    } else {
    graphics::plot(pscores, col = 0, pch = 16, ...,
        xlim = range(pscores[,1])*1.1, ylim =range(pscores[,2])*1.1,
        xlab = paste0("PC", dim[1], " ", percVar[dim[1]],"%"),
        ylab = paste0("PC", dim[2], " ", percVar[dim[2]],"%"))
    graphics::text(pscores, labels = rownames(pscores),
        col = pchCol, cex=cex.fac, xpd=FALSE)
    }

    if (!is.null(class) & legend != "none") {
        graphics::legend(legend, legend = levels(class), bg = "gray95",
            col = classCol, cex =.75, pch=16, text.col = classCol)
    }

    # add error labels
    if (!is.null(classConfusion) & !is.null(class)) {
        classConfusion <- classConfusion[keepN]
        highl <- class != classConfusion
        highl[is.na(highl) | is.na(class)] <- FALSE

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
    invisible(pemat)
}
