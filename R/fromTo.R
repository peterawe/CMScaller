#' translation between human gene identifiers
#' @description Translation between human gene identifiers. For input keys with
#' multiple valid matches, \bold{only first match is returned} (warning is
#' issued). If no \code{id.in} is provided, function will make an informed
#' guess.
#' @export
#' @param key a character vector, gene identifiers to be translated.
#' @param id.in a character, \code{colnames(anno.orgHs)} are valid options.
#' @param id.out a character, \code{colnames(anno.orgHs)} are valid options.
#' @param rough logical, whether to search for more than one match. If
#' @param all logical, whether to include all possible matches.
#' \code{length(key) > 10000}, setting \code{rough=TRUE} increases
#' speed by >100x, but also suppresses warnings.
#' @param verbose logical, whether console information are to returned.
#' @details Gene identifiers are from the \code{org.Hs.eg.db} package while
#' \emph{biotype} were retrieved from \code{biomaRt} using
#' \code{useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"}.
#' @return A character vector of same length as \emph{key} with gene
#' identifiers as specified.
#' @note Warnings are issued for keys with none or more than one hit. For keys
#' with more than one possible match, only first match is returned (thus
#' \emph{rough}).
#' @seealso \code{\link{anno.orgHs}}, \code{\link[org.Hs.eg.db]{org.Hs.eg.db}}
#' @references Marc Carlson (). org.Hs.eg.db: Genome wide annotation for Human.
#' R package version 3.2.3.
#' @references Steffen Durinck et al. (2009) Mapping identifiers for the
#' integration of genomic datasets with the R/Bioconductor package. Nature
#' Protocols 4, 1184-1191.
#' @examples
#'  fromTo(rownames(crcTCGAsubset)[1:50], "entrez", "symbol")
#'  head(anno.orgHs) # colnames indicate valid keys
#'  fromTo("AKT3", "symbol", "ensg") # expect one
#'  fromTo("AKT3", "symbol", "ensg", all = TRUE) # expect two
#'  # number of ids with multiple mappings
#'  sum(table(anno.orgHs$symbol) > 1)
#'  sum(table(anno.orgHs$entrez) > 1)
#'  sum(table(anno.orgHs$ensg) > 1)
#'  anno.orgHs[anno.orgHs$symbol == "AKT3",]
fromTo = function(key = NULL, id.in = NULL, id.out = "symbol",
                  verbose = getOption("verbose"),
                  rough = FALSE, all = FALSE) {

    #### checkInput ###########################################################

    if (sum(is.na(key)) > 0 & (rough == FALSE & all == FALSE))
        stop("NAs in key only accepted if rough=TRUE", call. = FALSE)

    # guess input id and use entrez as output in case of id.in == id.out
    if (is.null(id.in)) {
        mm <- apply(CMScaller::anno.orgHs, 2, function(x) sum(key %in% x))
        id.in <- colnames(CMScaller::anno.orgHs)[which.max(mm)]
        if (id.in == id.out) id.out <- "symbol"
        if (id.in == id.out) id.out <- "entrez"
        if (verbose == TRUE)
            message(paste0(id.in, " guessed as id.in, id.out is ", id.out))
    } else {
        if (!id.in %in% colnames(CMScaller::anno.orgHs)) stop("invalid id.in")
    }

    if (!id.out %in% colnames(CMScaller::anno.orgHs)) stop("invalid id.out")

    #### matchInform ##########################################################

    # but counts number of NAs and multi-matches
    if (all == TRUE & rough == FALSE) {
        tab <- CMScaller::anno.orgHs[!duplicated(paste0(
            CMScaller::anno.orgHs[,id.in], CMScaller::anno.orgHs[,id.out])) &
                CMScaller::anno.orgHs[,id.in] %in% key,c(id.in,id.out),drop = FALSE]

        # match ids and check for NAs and no-hits
        mm <- lapply(as.character(key), function(x) which(tab[,id.in] %in% x))
        res <- tab[unlist(mm), id.out]
        if ( length(res) != length(mm))
            message(paste0("length(key)=", length(key),
                           "; length(res)=", length(res)))

    }

    if (rough==TRUE & all == TRUE)  stop ("rough and all can not both be TRUE")

    if (rough==FALSE & all == FALSE) {

        # reduce search space -> matching is slow
        tab <- CMScaller::anno.orgHs[!duplicated(paste0(
            CMScaller::anno.orgHs[,id.in], CMScaller::anno.orgHs[,id.out])) &
            CMScaller::anno.orgHs[,id.in] %in% key,c(id.in,id.out),drop = FALSE]

        if (nrow(tab) == 0) stop("key not found - check id.in",call. = FALSE)

        # match ids and check for NAs and no-hits
        mm <- lapply(as.character(key), function(x) which(tab[,id.in] %in% x))
        isNA <- sum(sapply(mm, length) == 0)
        isMM <- sum(sapply(mm, length) > 1)

        # replace
        mm[sapply(mm, length) == 0] <- NA
        res <- tab[sapply(mm, "[[", 1), id.out]

        # return warnings
        if (isNA > 0)
            warning(paste0(isNA, " NA's returned"),
                           call. = FALSE)
        if (isMM > 0)
            warning(paste0(isMM, " multi-match; only one id/key returned!"),
                           call. = FALSE)
    }

    if (rough == TRUE & all == FALSE) {
    # matchFast ###############################################################
        res <- CMScaller::anno.orgHs[,id.out][
                match(key, CMScaller::anno.orgHs[,id.in])]
    }

    # sanity check
    if (length(key) != length(res) & all == FALSE)
        stop ("something is wrong!")

    # returnResults ###########################################################
    return(res)
}

