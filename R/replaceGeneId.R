#' replace matrix gene identifiers
#' @description Translate matrix rownames between gene identifiers.
#' @export
#' @param emat a numeric matrix with row features and sample columns.
#' @param id.in a character, gene symbol, entrez and ENSG are valid options.
#' @param id.out a character, gene symbol, entrez and ENSG are valid options.
#' @details Gene identifiers are from the
#' \code{\link[org.Hs.eg.db]{org.Hs.eg.db}} package.
#' @return A matrix with rownames gene identifier as specified in id.out.
#' Rows will be ordered by decreasing standard deviations.
#' \code{\link{make.unique}} are used to handle NA and duplicates in output.
#' @seealso \code{\link{anno.orgHs}}, \code{\link{fromTo}},
#' \code{\link[org.Hs.eg.db]{org.Hs.eg.db}}, \code{\link[stats]{sd}}
#' \code{\link{make.unique}}
#' @references Marc Carlson (). org.Hs.eg.db: Genome wide annotation for Human.
#' R package version 3.2.3.
#' @references Steffen Durinck et al. (2009) Mapping identifiers for the
#' integration of genomic datasets with the R/Bioconductor package. Nature
#' Protocols 4, 1184-1191.
#' @examples
#' # example data
#' x.entrez <- Biobase::exprs(crcTCGAsubset)
#' # rownames are entrez
#' head(rownames(x.entrez))
#' # translate to gene symbols
#' x.symbol <- replaceGeneId(x.entrez, id.in="entrez", id.out="symbol")
#' head(rownames(x.symbol))
#' # output matrix has fewer rows and row ordering is different
#' dim(x.entrez)
#' dim(x.symbol) # fewer rows
#' # x.entrez is not equal to x.entrez2
#' x.entrez2 <- replaceGeneId(x.symbol, id.in="symbol", id.out="entrez")
#' dim(x.entrez2)
replaceGeneId <- function(emat, id.in="symbol", id.out="entrez", verbose=TRUE) {
    if (class(emat) == "ExpressionSet") {
        emat <- suppressPackageStartupMessages(Biobase::exprs(emat))
    }
    if (class(emat) == "data.frame") emat <- as.matrix(emat)
    if (is.vector(emat)) emat <- matrix(emat, dimnames = list())
    if (is.null(rownames(emat))) stop("missing rownames(emat)")
    row.sd <- apply(emat, 1, stats::sd, na.rm=TRUE)
    emat <- emat[order(row.sd, decreasing = TRUE),]

    key.in  <- rownames(emat)
    key.out <- fromTo(key.in, id.in=id.in, id.out=id.out, rough=TRUE)

    rownames(emat) <- make.unique(key.out)

    return(emat)
}
