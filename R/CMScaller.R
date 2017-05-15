#' CMS classification
#' @description Cancer-cell intrinsic CMS classification based on pre-defined
#' subtype templates.
#' @export
#' @param emat a numeric expression matrix with sample columns, gene rows and
#' Entrez rownames. Microarray data should be normalized and log2-transformed.
#' For RNA-seq data, raw counts or RSEM values could be used directly by setting
#' \code{RNAseq=TRUE}.
#' @param RNAseq a logical, set to TRUE if emat is untransformed, non-normalized
#' sequencing counts or RSEM values.
#' @param nPerm an integer, number of permutations for \eqn{p}-value
#' estimation.
#' @param seed an integer, for \eqn{p}-value reproducibility.
#' @param pValue a p-value setting prediction confidence threshold.
#' @param verbose a logical, whether console messages are to be displayed.
#' @param doPlot a logical, whether to produce prediction \code{\link{subHeatmap}}.
#' @details \code{CMScaller} provides classification based on pre-defined
#' cancer-cell intrinsic CMS templates. The core algorithm is Nearest Template
#' Prediction (NTP) algorithm as proposed by Yujin Hoshida (2010). See
#' \code{\link{ntp}} for further details.
#' @note genes with missing values are discarded.
#' @return a data frame with class predictions, template distances,
#' \eqn{p}-values and false discovery rate adjusted \eqn{p}-values
#' (\code{\link[stats]{p.adjust}}). Rownames equal emat
#' colnames.
#' @seealso \code{\link{ntp}}, \code{\link{templates.CMS}},
#' \code{\link{corCosine}}
#' @references Hoshida, Y. (2010). Nearest Template Prediction: A
#' Single-Sample-Based Flexible Class Prediction with Confidence Assessment.
#' PLoS ONE 5, e15543.
#' @references Guinney J, Dienstmann R, Wang X, de Reynies A, Schlicker A,
#' Soneson C, et al. The consensus molecular subtypes of colorectal cancer.
#' Nat Med. 2015;21:1350-6.
#' @examples
#' res <- CMScaller(crcTCGAsubset, RNAseq=TRUE)
#' head(res)
#' hist(res$p.value)
CMScaller <- function(emat, RNAseq=FALSE, nPerm=1000, seed=42,
                    pValue=0.1, doPlot=TRUE, verbose=TRUE) {

    # checkInput ##############################################################

    # check datatype input and try to coerce to matrix
    if (class(emat) == "ExpressionSet") {
        emat <- suppressPackageStartupMessages(Biobase::exprs(emat))
    }
    if (class(emat) == "data.frame") emat <- as.matrix(emat)
    if (is.vector(emat)) emat <- matrix(emat, dimnames = list())
    if (is.null(rownames(emat))) stop("missing Entrez id rownames(emat)")

    if (ncol(emat) < 30) warnings("few samples - low prediction confidence")

    # quantile normalize RNA-seq data
    if (isTRUE(RNAseq)) {
        if (isTRUE(verbose)) message("performing quantile normalization...")
        emat <- log2(limma::normalizeBetweenArrays(emat)+.25) # + pseudo-count
    }

    # sanity check I - whether input data is log2-transformed
    emat.max <- abs(max(emat, na.rm = TRUE))
        if (emat.max > 25) {
            islog <- " <- check normalization and log-transformation"
            warning(paste0("emat max=", signif(emat.max,1), islog),
                    call. = FALSE)
        }

    # sanity check II - whether rownames appear to be Entrez ids
    mm <- sum(is.na(fromTo(rownames(emat), rough=TRUE)))/nrow(emat)
    if (mm > 0.15) warning ("check whether rownames(emat) are Entrez ids")

    # scale and center data, basically a wrapper for scale() function
    emat <- ematAdjust(emat)

    # ntpPredict ##############################################################

    res <- ntp(emat, CMScaller::templates.CMS, seed=seed,
                doPlot=doPlot, verbose=verbose)
    res <- subSetNA(res, pValue=pValue, verbose=verbose)

    # output ##################################################################

    # sanity check III - whether any FDR-values are above .1
    if (nPerm > 500) if (min(res$FDR) > .1)
        warning("low-confidence predictions - check input")

    return(res)
}
