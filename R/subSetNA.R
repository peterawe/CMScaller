#' sets NTP results to NA based on p-value cutoff
#' @export
#' @description sets \code{ntp} predictions to NA based on \eqn{p}-value or FDR
#' thresholds
#' @param res a data frame, result from \code{\link{ntp}} function.
#' @param pValue a numeric, predictions with higher p-values are set to NA.
#' @param FDR a numeric, predictions with higher FDR are set to NA.
#' @param verbose logical, whether console messages are to be displayed.
#' @details Replaces low-confidence predictions with NA.
#' @return A data.frame where res$prediction's are set to NA based on
#' \eqn{p}-value or FDR threshold.
#' @examples
#' emat <- ematAdjust(crcTCGAsubset, normMethod = "quantile")
#' res <- ntp(emat, templates.CMS, doPlot=TRUE, nPerm=100)
#' res <- subSetNA(res, pValue=.1)
#' res <- subSetNA(res, FDR=.1)
#' table(res$prediction, useNA="always")
subSetNA <- function(res, pValue = 1, FDR = 1, verbose = TRUE)
    {

    # checkInput #####
    if(is.null(res$prediction) | is.null(res$p.value) | is.null(res$FDR))
        { stop("error in res, check input") }

    # reset in case of prior resets
    classN <- levels(res$prediction)
    isDist <- grepl("d\\.", colnames(res))
    res$prediction <- classN[apply(res[,isDist], 1, which.min)]
    res$prediction <- factor(res$prediction, levels = classN)

    # process #####
    setNA <- res$p.value > pValue | res$FDR > FDR
    res$prediction[setNA] <- NA
    if (isTRUE(verbose))
        message(paste0(sum(setNA), "/", nrow(res), " samples set to NA"))
    # returnOutput #####
    return(res)
}
