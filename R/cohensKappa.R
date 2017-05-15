#' Cohen's kappa coefficient
#' @export
#' @description Calculates Cohen's \eqn{\kappa} (kappa) coefficient for two
#' equal length vectors.
#' @param x vector, coerced to factor.
#' @param y vector, coerced to factor.
#' @param adjusted logical. If \code{adjusted=FALSE}, the unadjusted accuracy
#' is returned instead.
#' @details Cohen's kappa coefficient measures the agreement between two
#' categorical vectors. Zero indicates that the agreement is no better than
#' expected by chance, while a value of one indicates perfect correspondence.
#' @note Cases with NA in either input vector are currently ignored.
#' @examples
#' x <- crcTCGAsubset$CMS
#' cohensKappa(x, x) # equals 1
#' replicate(10, cohensKappa(x, sample(x))) # expectation = 0
#' replicate(10, cohensKappa(x, sample(x), adjusted=FALSE)) # expectation > 0
cohensKappa <- function(x,y, adjusted=TRUE) {
    # input
    common.levels <- unique(c(as.character(x), as.character(y)))
    x <- factor(x, levels = common.levels)
    y <- factor(y, levels = common.levels)
    if (length(x) != length(y))
        stop ("x and y must be of equal length", call. = FALSE)
    # calculate
    tab <- base::table(x,y)
    observed.acc <- sum(diag(tab))/sum(tab)

    if (isTRUE(adjusted)) {
        marginal.freq.x <- rowSums(tab)/sum(tab)
        marginal.freq.y <- colSums(tab)/sum(tab)
        expected.acc <- sum(marginal.freq.x*marginal.freq.y)
        # Cohen's ĸ definitions
        return( 1-((1-observed.acc) / (1-expected.acc)) )
    } else {
        # simple accuracy
        return(observed.acc)
    }
}
