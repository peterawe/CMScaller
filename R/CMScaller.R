#' CMS classification
#' @description Cancer-cell intrinsic CMS classification based on pre-defined
#' subtype templates.
#' @export
#' @param emat a numeric expression matrix with sample columns, gene rows and
#' Entrez rownames. Microarray data should be normalized. For RNA-seq data,
#' counts or RSEM values could be used directly by setting \code{RNAseq=TRUE}.
#' @param templates a data frame with two columns; \emph{class} (coerced to
#' factor) and \emph{probe} (coerced to character).
#' @param rowNames a character, either "entrez" (NCBI Entrez),
#' "symbol" (HGNC symbol) or "ensg" (Ensembl). If set to other than "ensg",
#' \code{\link{replaceGeneId}} is used to translate \code{rownames(emat)}.
#' @param RNAseq a logical, set to TRUE if emat is untransformed, non-normalized
#' sequencing counts or RSEM values.
#' @param nPerm an integer, number of permutations for \eqn{p}-value
#' estimation.
#' @param seed an integer, for \eqn{p}-value reproducibility.
#' @param FDR a false discovery rate, sets prediction confidence threshold.
#' @param verbose a logical, whether console messages are to be displayed.
#' @param doPlot a logical, whether to produce prediction \code{\link{subHeatmap}}.
#' @details \code{CMScaller} provides classification based on pre-defined
#' cancer-cell intrinsic CMS templates. If \code{RNA-seq=TRUE}, a pseudocount
#' of 0.25 is added, matrix log2-transformed and quantile normalized
#' (\code{\link[limma]{normalizeQuantiles}}) prior to scaling/centering and
#' prediction. The core algorithm is the Nearest Template Prediction (NTP)
#' algorithm as proposed by Yujin Hoshida (2010). See \code{\link{ntp}} for
#' further details.
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
CMScaller <- function(emat, templates=CMScaller::templates.CMS,
                    rowNames="entrez",
                    RNAseq=FALSE, nPerm=1000, seed=NULL,
                    FDR=0.05, doPlot=TRUE, verbose=TRUE) {

    # checkInput ##############################################################

    # check datatype input and try to coerce to matrix
    if (class(emat)[1] == "ExpressionSet") {
        emat <- suppressPackageStartupMessages(Biobase::exprs(emat))
    }
    if (class(emat)[1] == "data.frame") emat <- as.matrix(emat)
    if (is.vector(emat)) emat <- matrix(emat, dimnames = list())
    if (is.null(rownames(emat))) stop("missing Ensembl id rownames(emat)")

    if (ncol(emat) < 30) warnings("few samples - high prediction variance",
                                call.=FALSE)

    if (rowNames != "entrez") {
        if (!rowNames %in% c("symbol", "ensg"))
            stop("invalid rowNames, must be either entrez, symbol or ensg")
            emat <- replaceGeneId(emat, id.in=rowNames, id.out="entrez")
    }

    # log2-transform and quantile normalize RNA-seq data
    if (isTRUE(RNAseq)) {
        if (isTRUE(verbose))
            message("performing log2-transform and quantile normalization...")
        emat <- limma::normalizeQuantiles(log2(emat+.25))
    }

    # sanity check - whether rownames appear to be Entrez ids
    is.na.rows <- is.na(fromTo(rownames(emat), rough=TRUE))
    mm <- sum(is.na.rows)/nrow(emat)
    if (mm > 0.15) {
        message (paste0(sum(is.na.rows),"/",nrow(emat),
                    " rownames(emat) failed to match to human gene identifiers"))
        warning (paste0("verify that rownames(emat) are ", rowNames),
                call.=FALSE)
    }

    # scale and center data, basically a wrapper for scale() function
    emat <- ematAdjust(emat)

    # ntpPredict ##############################################################

    res <- ntp(emat, templates, seed=seed, nPerm=nPerm,
                doPlot=doPlot, verbose=verbose)
    res <- subSetNA(res, FDR=FDR, verbose=verbose)

    # output ##################################################################

    # sanity check III - whether any FDR-values are above .1
    if (nPerm > 500) if (min(res$FDR) > .1)
        warning("low-confidence predictions - check input",call.=FALSE)

    return(res)
}
