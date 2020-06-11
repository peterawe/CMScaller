#' replace matrix gene identifiers
#' @description Translate matrix rownames between gene identifiers.
#' @export
#' @param emat a numeric matrix with row features and sample columns.
#' @param id.in a character, gene symbol, entrez and ENSG are valid options.
#' @param id.out a character, gene symbol, entrez and ENSG are valid options.
#' @details Gene identifiers are from the
#' \code{\link[org.Hs.eg.db]{org.Hs.eg.db}} package.
#' @return A matrix with rownames gene identifier as specified in id.out.
#' For duplicates, rownames for rows with lower standard deviation are added a
#' number to make unique.
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
#' library(Biobase)
#' x.entrez <- Biobase::exprs(crcTCGAsubset)
#' # rownames are entrez
#' head(rownames(x.entrez))
#' # translate to gene symbols
#' x.symbol <- replaceGeneId(x.entrez, id.in="entrez", id.out="symbol")
#' head(rownames(x.symbol))
#' x.entrez2 <- replaceGeneId(x.symbol, id.in="symbol", id.out="entrez")
#' # translations are not cycle consistent
#' table(rownames(x.entrez2) == rownames(x.entrez))
#' # matrix values are not changed
#' all(x.entrez == x.entrez2)
replaceGeneId <- function(emat, id.in="symbol", id.out="entrez") {
    if (class(emat)[1] == "ExpressionSet") {
        emat <- suppressPackageStartupMessages(Biobase::exprs(emat))
    }
    if (class(emat)[1] == "data.frame") emat <- as.matrix(emat)
    if (is.vector(emat)) emat <- matrix(emat, dimnames = list())
    if (is.null(rownames(emat))) stop("missing rownames(emat)")
    row.sd <- apply(emat, 1, stats::sd, na.rm=TRUE)

    key.in  <- rownames(emat)
    key.out <- fromTo(key.in, id.in=id.in, id.out=id.out, rough=TRUE)
    key.na <- is.na(key.out)
        message (paste0(sum(key.na), "/", length(key.na),
                        " rownames [NA.number] (no valid translation)"))
        message (paste0(sum(duplicated(stats::na.omit(key.out))), "/", length(key.na),
                        " rownames [id.number] (translation gives duplicates)"))

    # ordering by row.sd -> most informative probe are not padded with number
    key.out[is.na(key.out)] <- "NA"
    key.out <- make.unique(key.out[order(row.sd, decreasing = TRUE)])
    key.out <- key.out[match(key.in, key.in[order(row.sd, decreasing = TRUE)])]
    key.out <- make.unique(key.out)
    rownames(emat) <- key.out
    return(emat)
}

