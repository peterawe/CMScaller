#' nearest template prediction templates
#' @export
#' @description Reformats list of gene sets or results from
#' \code{\link{subDEG}} into prediction templates for \code{\link{ntp}}.
#' @param geneSets list of vectors with gene identifiers (coerced to
#' character).
#' @param resDEG logical, if \code{geneSets} is result from \link{subDEG}.
#' @param lfc numeric, log2fold-change threshold (only applicable if
#' \code{resDEG=TRUE}).
#' @param padj numeric, \eqn{p_{adj}}-value threshold (only applicable if
#' \code{resDEG=TRUE}).
#' @param topN integer, (maximum) number of genes per class.
#' @param verbose logical, whether console messages are to be displayed.
#' @seealso \code{\link{ntp}}, \code{\link{subDEG}}
#' @return
#' A data.frame formated as templates for \code{\link{ntp}}.
#' @examples
#' deg <- subDEG(crcTCGAsubset, crcTCGAsubset$CMS, doVoom=TRUE, sortBy="B")
#' tmp <- ntpMakeTemplates(deg, topN=25)
#' table(tmp$class)
ntpMakeTemplates <- function(geneSets, resDEG=TRUE,
                            lfc = 1, padj=.1, topN=NULL,
                            verbose=getOption("verbose", FALSE)) {

    if (isTRUE(resDEG)) {
        if (class(geneSets) == "list") {
            geneSets <- lapply(geneSets, function(deg) {
                rownames(deg)[which(deg$logFC > lfc & deg$adj.P.Val < padj)]
            })

        if (!is.null(topN)) geneSets <- lapply(geneSets, utils::head, n=topN)

        } else {
            subtypes <- attr(geneSets, "contrast")
            geneSets <- split(geneSets, geneSets$logFC < 0)

            if (!is.null(topN)) {
                geneSets <- lapply(geneSets, utils::head, n=topN)
                geneSets <- lapply(geneSets, rownames)
            } else {
                geneSets <- lapply(geneSets, function(deg) {
                    rownames(deg)[abs(deg$logFC) > lfc & deg$adj.P.Val < padj]
                })
            }
            names(geneSets) <- subtypes
        }
    }

    if (!all(unlist(lapply(geneSets, typeof)) %in%
             c("factor", "character", "numeric"))) stop ("error input")
    if (is.null(names(geneSets))) stop ("no genesets names")
    if (length(geneSets) < 2) stop ("need at least two genesets")
    if (min(sapply(geneSets,length)<=5)) warning ("few features per class")

    class <- make.names(names(geneSets))
    probe <- unlist(lapply(geneSets, as.character))
    ngenes <- sapply(geneSets,length)
    templates <- data.frame(probe = unlist(probe),
                            class = factor(rep(class, ngenes), levels=class),
                            stringsAsFactors = FALSE,
                            check.names = FALSE)

    feat.class <- paste(range(table(templates$class)),collapse = "-")

    if (isTRUE(verbose)) message(paste0(feat.class, " features/class\n",
            "classes: ", paste(class, collapse=" ")))
    rownames(templates) <- NULL
    return(templates)
}
