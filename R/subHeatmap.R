#' ntp heatmap
#' @description Heatmap with marker genes and class labels.
#' @export
#' @param emat a numeric matrix with row features and sample columns.
#' @param res a data frame, result from \code{ntp} function.
#' @param templates data frame with two columns; class (factor or numeric) and
#' probe (character).
#' @param keepN a logical or numeric vector specifying which samples to keep.
#' @param classCol a character vector specifying class colors.
#' @param labRow a character specifying \code{templates} column used for labeling
#' features.
#' @param heatCol a character vector specifying heatmap colors (must be an even
#' number).
#' @param ... additional arguments passed to \code{\link[graphics]{image}}.
#' @details \code{subHeatmap} provides a simple heatmap with samples and genes
#' ordered according to subtype. The column bar indicates prediction
#' \eqn{p}-values.
#' @return a raster heatmap to visualize template feature expression.
subHeatmap <- function(emat, res, templates, keepN = TRUE,
                       labRow = NULL,
                       classCol = getOption("subClassCol"),
                       heatCol = NULL,...) {

    # checkInput ##############################################################

    if (!all(colnames(emat) == rownames(res))) stop("check if emat/res match")

    # prepareInput ############################################################

    emat <- emat[,keepN];res <- res[keepN,]
    class <- res$prediction
    K <- nlevels(class)

    if (!is.null(labRow )) {
        old.par <- graphics::par()
        on.exit(graphics::par(mar=old.par$mar))
        graphics::par(mar=old.par$mar+c(0,0,0,3))
    }

    # prepare templates
    templates <- templates[!duplicated(templates$probe),]
    rownames(templates) <- templates$probe

    # prepare res
    isDist <- grepl("^d\\.", colnames(res))
    res$prediction <- factor(levels(res$prediction)
                            [apply(res[,isDist], 1, which.min)])

    if (K == 2) {
        isClass2 <- res$prediction == levels(class)[2]
        res$p.value[isClass2] <- -res$p.value[isClass2]
    }

    # match emat and templates and order colums after prediction and p.value
    N <- order(res$prediction, res$p.value)
    P <- intersect(rownames(emat), rownames(templates))

    templates <- templates[P,]
    templates <- templates[order(templates$class),]
    emat <- ematAdjust(emat[rownames(templates),N])

    # colorBars ###############################################################

    if (is.null(heatCol)) heatCol <- subData[["hmCol"]]

    # prettyBreaks ############################################################
    pMax = 3
    breaksLeft = seq(-pMax,-.25, length = length(heatCol)/2)
    breaks <- c(breaksLeft, 0, abs(rev(breaksLeft)))

    # necessary to avoid color overflow
    emat[emat > pMax] <- pMax; emat[emat < -pMax] <- -pMax

    # makeHeatmap #############################################################
    xx <- seq(0, 1, length.out=ncol(emat)+1)
    yy <- seq(0, 1, length.out=nrow(emat)+1)

    graphics::image(x = xx, y = yy, z = t(emat),
                yaxt="n", xaxt="n",useRaster=TRUE,
                col = heatCol, breaks = breaks,
                xlab = "class predictions", ylab = "template features", ...)

    # make sample bar
    xb <- cumsum(c(0,(sapply(split(res$prediction[N], res$prediction[N]),
                             length)/length(N))))
    xl <- xb[-(K+1)]
    xr <- xb[-1]

    yy <- line2user(line=1:0, side=1)
    yy <- seq(yy[1], yy[2], length=3)

    graphics::rect(xleft = xl, xright = xr, ybottom = yy[1], ytop=yy[2],
                   col=classCol, xpd=TRUE, border=FALSE, lwd = .75)

    xxl <- xl+(xr-xl)/2
    graphics::text(xxl, yy[1], pos=1,levels(class), adj=0,xpd=TRUE, cex=.75)

    # make p-value trace if applicable
    if (min(res$p.value) < 1) {
        py <- abs(res$p.value[N])
        px <- (seq_along(N)/(length(N)))-1/length(N)/2
        py <- (py*(yy[2]-yy[1]))+yy[1]
        graphics::rect(xleft = xx[-length(xx)], xright = xx[-1],
                       ybottom=yy[1], ytop=py, border = FALSE,
                       xpd=NA, col = "white", lwd=2)
        graphics::rect(xleft = xl, xright = xr, ybottom = yy[1], ytop=yy[2],
                       col=NA, xpd=TRUE, lwd = .75)

        graphics::text(0,yy[1], pos=2, expression(italic(p)-value), xpd=TRUE, cex=.75)
        graphics::text(x=1,y=c(yy[1],yy[2]), c(0,1), pos=c(4),xpd=TRUE,cex=.75)
        graphics::text(4,2, paste(bquote(italic(p)-value)))
    }

    # make feature bar
    xx <- line2user(line=1:0, side=2)
    xx <- seq(xx[1], xx[2], length=3)

    yy <- c(0, cumsum(sapply(split(templates$probe, templates$class),
                            length)/length(templates$probe)))
    yb <- yy[-(K+1)]
    yt <- yy[-1]
    graphics::rect(xleft=xx[1],xright=xx[2],ybottom=yb,ytop=yt,
                   col=classCol, xpd=NA, lwd = .75)

    if (!is.null(labRow )) {
        textfun <- function(..., cex.text)
            graphics::text(..., xpd=TRUE, cex=.75)
        id <- templates[,labRow]
        yy <- seq(0+(1/nrow(templates)/2), 1-(1/nrow(templates)/2),
                  length.out=nrow(templates))
        textfun(1, yy, id, pos=4, adj=1)
    }
}
