#' MA plot
#' @export
#' @description Mean versus average (MA) plot of output from \link{subDEG}.
#' @param deg output from \link{subDEG} which may be either a list or
#' data.frame.
#' @param geneID character, column names to use for labeling top hits.
#' @param lfc numeric absolute \eqn{log2(fold-change)} threshold.
#' @param padj numeric, adjusted \eqn{p}-value threshold.
#' @param ave numeric, average expression threshold.
#' \code{length(class)==ncol(emat)}.
#' @param topN an integer, number of genes to label in plot (selected by
#' largest absolute \eqn{fold-change}).
#' @param classCol a character vector specifying class colors.
#' @param cexText numeric, scaling factor for labels for top hits.
#' @param ... additional arguments passed to \code{\link[graphics]{plot}} or
#' \code{\link[graphics]{smoothScatter}} or if \pkg{KernSmooth} is available
#' (currently only main, xlab and ylab used).
#' @return a MA-plot
#' @examples
#' deg <- subDEG(crcTCGAsubset, crcTCGAsubset$CMS, doVoom=TRUE, sortBy="B")
#' subMA(deg, geneID="symbol", lfc=3)
subMA <- function(deg, geneID = "rownames",
                        lfc = log2(2), padj = .01, ave = 1,
                        topN = 10, cexText=1,
                        classCol = getOption("subClassCol"),...) {

    ddd <- list(...)
    if (!is.null(ddd$main)) main <- ddd$main else main=""
    if (!is.null(ddd$xlab)) xlab <- ddd$xlab else
        xlab=expression(mean(log[2](signal)))
    if (!is.null(ddd$ylab)) ylab <- ddd$ylab else
        ylab=expression(-log[2](foldchange))

    plotVolc <- function(deg, clCol = classCol) {
        y <- deg$logFC
        x <- deg$AveExp
        p <- deg$adj.P.Val


    if (geneID=="rownames") geneID <- rownames(deg) else geneID <- deg[,geneID]

        # plot background
        xRange <- c(min(x),  max(x, na.rm=TRUE)*1.05)
        yRange <- c(-max(y, na.rm=TRUE)*1.25, max(y, na.rm=TRUE)*1.2)

        if (length(x) < 3e3) {
            graphics::plot(x, y, xlim=xRange, col="gray", cex=.5,
                            main=main, xlab=xlab, ylab=ylab)
        } else {
            if (!packageExists("KernSmooth")) {
                graphics::plot(x, y, xlim=xRange, col="gray", cex=.5,
                               main=main, xlab=xlab, ylab=ylab)
            } else {
                #graphics::smoothScatter(x, y, xlim=xRange, nrpoints=0,
                #            main=main, xlab=xlab, ylab=ylab)
            }
        }
        graphics::abline(h=0)

        # add features crossing threshollds
        ff <- which(p < padj & abs(y) > lfc)
        if (sum(ff) > 1) {
            cc <- ifelse(y > 0, clCol[1], clCol[2])
            graphics::points(x[ff], y[ff], col=cc[ff], cex=.5)
            }

            # label top-DEGs
            if (topN > 0) {
                ff2 <- seq_len(nrow(deg)) %in% ff & abs(deg$AveExpr) > ave
                dn <- rank(y[ff2]) <= topN
                up <- rank(-y[ff2]) <= topN
                top <- dn | up
                if (sum(top) > 0)
                    graphics::text(x[ff2][top], y[ff2][top],geneID[ff2][top], cex=cexText)
            }

    return(data.frame(x, y, geneID,
                        row.names = rownames(deg), stringsAsFactors = FALSE))
    }

    mtextFun <- function(...) graphics::mtext(...,side=3, cex=.67, line=0)

    if(class(deg) == "list") {
        K <- length(deg)
        snk <- lapply(seq_len(K), function(k) {
            plotVolc(deg[[k]],
                     clCol = c(classCol[k %% length(classCol)], "gray"))
            #mtextFun(paste(names(deg)[k],"up "), adj=1)
            #mtextFun(paste(" ",names(deg)[k],"down"), adj=0)
            })
    } else {
        snk <- plotVolc(deg)
        #mtextFun(paste(" ",attr(deg, "contrast")[1]), adj=1)
        #mtextFun(paste(attr(deg, "contrast")[2], " "), adj=0)
    }
    invisible(snk)
}
