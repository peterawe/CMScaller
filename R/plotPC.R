#' 2d plot of PCA variable loadings
#' @description Plots the loadings for the highest contributing variables (by euclidean distance).
#' @export
#' @param x list, PCA output (see example code)
#' @param dim a numeric vector of length 2, specifying which principal
#' components to plot.
#' @param n numeric, number of variables to include
#' @param entrez logical, if TRUE, translates variables from Entrez to Gene symbosl
#' @param ... additional arguments passed to \code{\link[graphics]{plot}}
#' @examples
#' library(CMScaller)
#' p <- subPCA(crcTCGAsubset[,1:10])
#' plotPC(p, entrez=TRUE)
#' @seealso \code{\link[stats]{prcomp}}, \code{\link{subPCA}}
plotPC <- function(x, dim=c(1:2), n=10, entrez=FALSE, ...) {
    offset <- 1.05
    # scale components
    x <- (x$rotation %*% diag(x$sdev))[,dim,drop=FALSE]
    # include only top-n highest euclidean length (radius)
    xsum <- apply(x^2, 1, sum)
    o <- utils::head(order(xsum, decreasing = TRUE), n)
    x <- x[o,,drop=FALSE]
    colnames(x) <- paste0("PC", dim)
    xmax=max(abs(x))*offset
    # prepare window
    graphics::plot(0,0,xlim=c(-xmax,xmax), ylim=c(-xmax,xmax),
                   xlab=paste0("PC", dim[1]), ylab=paste0("PC", dim[2]),
                   col=0, frame=FALSE, asp = 1, ...)
    plotrix::draw.circle(0,0,max(sqrt(xsum)), lty=2, border = "gray")

    if(entrez==TRUE) rownames(x) <- make.unique(CMScaller::fromTo(rownames(x)))
    # draw and label
    graphics::arrows(x0 = 0, y0 = 0, x[,1], x[,2], col="gray", length = .1)

    graphics::text(x[,1]*offset, x[,2]*offset, rownames(x), cex=.75,xpd=TRUE)
    invisible(x)
}

