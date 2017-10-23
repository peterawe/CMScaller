#' scatterplots of sample-templates distances
#' @export
#' @description Pairwise sample-templates distance plots with samples colored
#' according to predicted class.
#' @param res a data frame, result from \code{\link{ntp}} function.
#' @param classCol a character vector specifying class colors
#' @return a panel of plots of pair-wise template-distances.
#' @details Upper-panels show the Pearson correlations
#' between sample-templates distances. Diagonal panels show class and histogram
#' of distances. The vertical, dashed lines are the distance between two
#' uncorrelated vectors.
#' @note smaller distances indicate more similar samples.
#' @seealso \code{\link{ntp}}, \code{\link{corCosine}}, \code{\link{pairs}}
#' @examples
#'  res <- CMScaller(crcTCGAsubset, nPerm=1000, RNAseq=TRUE)
#'  subPairs(res)
subPairs <- function(res, classCol = getOption("subClassCol"))
    {


    ###########################################################################
    # process and output
    ###########################################################################

    p2hex.alpha <-  function(p) {
        sprintf("%02X", floor((1-p^2) * 255))
    }

    gray <- "#bbbbbb"
    pchCol <- classCol[res$prediction]
    pchCol[is.na(pchCol)] <- gray
    pchCol <- paste0(pchCol,p2hex.alpha(res$p.value))

    panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
    {
        usr <- graphics::par("usr"); on.exit(graphics::par(usr))
        graphics::par(usr = c(0, 1, 0, 1))
        r <- stats::cor(x, y)
        txt <- format(c(r, 0.123456789), digits = digits)[1]
        graphics::text(0.5, 0.5, txt)
    }

    panel.hist <- function(x, ...) {
        usr <- graphics::par("usr"); on.exit(graphics::par(usr))
        graphics::par(usr = c(usr[1:2], 0, 1.5) )
        h <- graphics::hist(x, plot = FALSE)
        breaks <- h$breaks; nB <- length(breaks)
        y <- h$counts; y <- y/max(y)
        graphics::rect(breaks[-nB], 0, breaks[-1], y)
        graphics::abline(v=sqrt(1/2 * (1 - (0))),
                        lty=2, xpd=FALSE,col=gray)
    }

    isDist <- grepl("^d\\.", colnames(res))
    maxDist <- max(res[,isDist],na.rm=TRUE)
    minDist <- min(res[,isDist],na.rm=TRUE)
    graphics::pairs(res[,isDist],
        pch = 16, col = pchCol,
        diag.panel = panel.hist,
        upper.panel = panel.cor,
        xlim = c(minDist,maxDist), ylim = c(minDist,maxDist))
    graphics::mtext("sample-template correlation distance", side = 1, line = 4, cex=.75)
    graphics::mtext("sample-template correlation distance", side = 2, line = 3, cex=.75)
}
