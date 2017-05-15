#' boxplot with data points
#' @description Boxplot with data points and labeled outliers.
#' @export
#' @param y a numeric vector.
#' @param x a factor, \code{length(x)==length(y)}.
#' @param labels a character vector, if NULL, \code{names(y)} is used.
#' @param classCol a character vector specifying class colors.
#' @param keepN logical or numeric indicating which samples to include.
#' @param points logical indicating whether points are to be added to boxes.
#' @param width numeric to adjust point stacking - increase to reduce width.
#' @param showOutliers logical indicating whether outliers are highlighted.
#' @param ... additional arguments to be passed to \code{graphics::boxplot}.
#' @note For data points, \eqn{y}-values are grouped and averaged - not exact!
#' @examples
#' isCMS4 <- rownames(crcTCGAsubset) %in%
#'     templates.CMS$probe[templates.CMS$class == 'CMS4']
#' y <- colSums(Biobase::exprs(crcTCGAsubset)[isCMS4,])
#' x <- crcTCGAsubset$CMS
#' subBoxplot(y, x, ylab = expression(sum(log2(RSEM))), main = 'CMS4 genes',
#'      keepN = which(x != "CMS1"), notch = TRUE)
subBoxplot = function(y, x, labels = NULL, keepN = TRUE,
                    classCol = getOption("subClassCol"),
                    points = FALSE, width = 20, showOutliers = FALSE,...)
    {

    if (is.null(labels)) {
            if (is.null(names(y))) labels <- seq_along(y)
                else labels <- names(y)
    }

    # modified from http://www.r-bloggers.com/labeled-outliers-in-r-boxplot/
    data <- data.frame(x,y=as.vector(y), labels)[keepN,,drop = FALSE]
    data <- droplevels(data)

    # quick fix to include NA in labels
    if (anyNA(data$x)) {
        if (is.factor(data$x)) {
            class.names <- levels(data$x)
        } else {
            class.names <- unique(data$x)
        }
        data$x <- as.character(data$x)
        data$x[is.na(data$x)] <- " NA"
        data$x <- factor(data$x, levels=c(class.names, " NA"))
    }

    # pch = "" is a hack to avoid double plotting!
    boxdata <- with(data,graphics::boxplot(data$y ~ data$x,
                                plot = TRUE, pch=ifelse(isTRUE(points),"",1),
                                col = classCol, ...))

    # add data points to boxplot (christmas)
    if (points == TRUE) {
        yList <- split(data$y, data$x)
            sapply(seq_along(yList), function(k) {
                yc <- sort(cut(yList[[k]], 50))
                # approximation!!! for pretty stacking
                yApprox <- sapply(strsplit(gsub("\\(|\\]", "", yc), "\\,"),
                        function(x) mean(as.numeric(x)))
                xt <- table(yApprox)
                xx <- sapply(seq_along(xt), function(i) {
                    x <- seq(0,length = xt[i])/ width
                    x-mean(x)+k
            })

        graphics::points(unlist(xx), yApprox, pch=16, col="lightgray")
        })

        # add text to boxplot
        if (length(boxdata$out) > 0 & showOutliers == TRUE) {
            for (i in seq_along(boxdata$group)) {
                graphics::text(boxdata$group[i], boxdata$out[i],
                    data$label[which(
                        data$y==boxdata$out[i] &
                        as.numeric(addNA(data$x,ifany=TRUE))==boxdata$group[i])],
                cex = 0.75)
        }
    }
    }
}
