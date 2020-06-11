#' limma/voom differential gene expression analysis
#' @export
#' @description Wrapper function to perform \code{\link[limma]{eBayes}} differential
#' gene expression analysis (DEG).
#' @param emat numeric matrix with row features and sample columns. For
#' \code{\link[Biobase]{ExpressionSet}}
#' \code{\link[Biobase]{exprs}} matrix is used and featureData used for
#' annotating results.
#' @param class a factor vector specifying subtypes compared.
#' \code{length(class)==ncol(emat)}.
#' @param batch a factor vector specifying additional level in design matrix.
#' @param keepN a logical or numeric vector specifying which samples to keep
#' (defaults to all).
#' @param feat data frame with feature data
#' @param doPairwise a logical, whether pairwise comparison are to be
#' performed (as opposed to class k against K-k).
#' @param doPlot a logical, not implemented.
#' @param returnTreat logical, if true, parameters lfc and padj are passed to \code{\link[limma]{topTreat}}, else \code{\link[limma]{eBayes}} are used and results sorted as indicated by sortBy parameter.
#' @param lfc numeric, only applicable if returnTreat=TRUE, if FALSE, all genes are returned
#' @param padj numeric, only applicable if returnTreat=TRUE, if FALSE, all genes are returned
#' @param doVoom a logical, indicating whether emat is sequencing count data.
#' @param normMethod a character, only used if doVoom=TRUE and passed to
#' @param ... additional arguments passed to \code{\link[edgeR]{calcNormFactors}}
#' \code{\link[edgeR]{calcNormFactors}} if element in c("TMM","RLE",
#' "upperquartile","none") or \code{\link[limma]{voom}} ("scale", "quantile",
#' "cyclicloess").
#' @param sortBy character passed to \code{\link[limma]{topTable}}
#' \code{sort.by} parameter.
#' @return a data frame (for two classes) or list of data frames (in case of
#' more than two classes) with results.
#' @references Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. Nucl. Acids Res. 2012;gks461. \url{http://nar.oxfordjournals.org/content/early/2012/05/24/nar.gks461}
#' @examples
#' deg <- subDEG(crcTCGAsubset, crcTCGAsubset$CMS, doVoom=TRUE, sortBy="P")
#' lapply(deg, head)
subDEG <- function(emat, class, batch=NULL, keepN=TRUE, doPairwise = FALSE,
                    doVoom = FALSE, normMethod = "quantile",
                    returnTreat=FALSE, lfc=log2(1.5), padj=0.05,
                    feat=NULL, doPlot=FALSE, sortBy="P", ...)

{

    if (!packageExists("limma"))
        stop (paste0("Function requires limma package available from:\n",
                     "http://www.bioconductor.org/packages/release/bioc/html/limma.html"))

    ###########################################################################
    #  test/fix/clean input
    ###########################################################################

    if (class(emat)[1] == "ExpressionSet") {
        if (is.null(feat)) feat <- Biobase::fData(emat)
        emat <- Biobase::exprs(emat)
    }

    # drop samples with unknown class label
    if (sum(is.na(class[keepN])) > 0 ) {
        message(paste0(sum(is.na(class[keepN])),
                       " samples with class or batch NA's excluded"))
        if (is.numeric(keepN)) keepN <- seq_len(ncol(emat)) %in% keepN
        keepN <- keepN & !is.na(class)
    }

    # coerce synthatically valid class labels
    class.levels <- levels(factor(class))
    class <- factor(make.names(class))
    if (!is.null(batch)) batch <- factor(make.names(batch))

    if (all(class.levels %in% levels(class)))
        class <- factor(class, levels = class.levels)

    # drop samples with unknown class label
    if (sum(is.na(class[keepN])) > 0 ) {
        message(paste0(sum(is.na(class[keepN])),
                       " samples with class or batch NAs excluded"))
        if (is.numeric(keepN)) keepN <- seq_len(ncol(emat)) %in% keepN
        if (!is.null(batch)) {
            keepN <- keepN & !is.na(class) & !is.na(batch)
        } else {
            keepN <- keepN & !is.na(class)
        }
    }

    requireNamespace("limma", quietly = TRUE)

    ###########################################################################
    #  BATCH == NULL: specify DEG modeling
    ###########################################################################

    if(is.null(batch)) {
        class <- droplevels(class[keepN])
        design <- stats::model.matrix(~0+class, data = class)
        colnames(design) <- levels(class)
        emat <- emat[,keepN]

        if (doVoom == TRUE) emat <- voomTransform(emat, normMethod, design, ...)

        fit <- limma::eBayes(limma::lmFit(emat, design))
        fit2 <- limma::eBayes(fit)

        # make comparisons
        mm <- levels(factor(class))
        if(doPairwise == TRUE) {
            mm = t(utils::combn(mm, 2))
            cont.gen = vector("list", length = nrow(mm))
            names(cont.gen) = apply(mm, 1, paste, collapse = "-")
            loopLength = nrow(mm)

        } else {
            cont.gen = vector("list", length = length(mm))
            names(cont.gen) = as.character(mm)
            loopLength = length(mm)
            if(loopLength == 2) loopLength = loopLength -1
        }

        #######################################################################
        #  DEG modeling
        #######################################################################

        for(m in seq_len(loopLength)) {
            if (doPairwise == TRUE) comparison =
                            paste(mm[m,1],mm[m,2], sep ="-")
            if (doPairwise == FALSE) comparison = paste0(mm[m], "-(",
                            paste(mm[-m], collapse= "+"),")/", length(mm)-1)
            contrasts2 <- limma::makeContrasts(contrasts = comparison,
                                               levels=design)
            contrFit2 <- limma::eBayes(limma::contrasts.fit(fit2, contrasts2))

            if (returnTreat == TRUE) {
                treatFit = limma::treat(limma::contrasts.fit(fit2, contrasts2),
                                        lfc = lfc)
                cont.gen[[m]] = limma::topTreat((treatFit),
                        number = Inf, p.value = padj, genelist = feat,
                        sort.by = sortBy)
            } else cont.gen[[m]] = limma::topTable(contrFit2, genelist = feat,
                        number = Inf, sort.by= sortBy)
            #comparison <- gsub("-.*", " vs rest", comparison)
            # if(isTRUE(doPlot))
            #     subVolcano(contrFit2, lfc = lfc, padj = padj)
        }


    ###########################################################################
    #  BATCH != NULL: specify DEG modeling
    ###########################################################################

    } else {
        class <- droplevels(class[keepN])
        batch <- batch[keepN]
        df <- droplevels(data.frame(group = class, batch = batch))
        design <- stats::model.matrix(~0+group+batch, data = df)
        colnames(design) <- gsub("group", "", colnames(design))
        emat <- emat[,keepN]

        if (doVoom == TRUE) emat <- voomTransform(emat, normMethod, design, ...)

        fit <- limma::lmFit(emat, design)
        fit2 <- limma::eBayes(fit)

        # make comparisons
        mm <- levels(df$group)
        loopLength <- length(mm)
        cont.gen <- vector("list", length = length(mm))
        names(cont.gen) <- mm

        if(doPairwise == TRUE) {
            mm = t(utils::combn(mm, 2))
            cont.gen = vector("list", length = nrow(mm))
            names(cont.gen) = apply(mm, 1, paste, collapse = "-")
            loopLength = nrow(mm)
        }

        #######################################################################
        #  DEG modeling
        #######################################################################

        for(m in seq_len(loopLength)) {
            if(doPairwise == TRUE) {
                contrasts1 = paste(mm[m,1],mm[m,2], sep ="-")
                contrasts2 = limma::makeContrasts(contrasts = contrasts1,
                                        levels = colnames(design))}
            if(doPairwise == FALSE) {
                contrasts1 = paste0(mm[m],"-(", paste0(mm[-m],
                                        collapse = "+"),")/",length(mm)-1)
                contrasts2 = limma::makeContrasts(contrasts =
                                        contrasts1, levels = colnames(design))
            }
            cont.title <- gsub("\\+.*", "...", contrasts1)
            contrFit2 <- limma::eBayes(limma::contrasts.fit(fit,
                                        contrasts = contrasts2))
            if(returnTreat == TRUE) {
                treatFit = limma::treat(limma::contrasts.fit(fit2,
                                                     contrasts2), lfc = lfc)
                cont.gen[[m]] = limma::topTreat((treatFit),
                                            number = Inf, p.value = padj,
                                            genelist = feat, sort.by = sortBy)
            } else cont.gen[[m]] = limma::topTable(contrFit2, genelist = feat,
                                            number = Inf, sort.by= sortBy)
            #comparison = gsub("-.*", " vs rest", comparison)
            # if(isTRUE(doPlot))
            #     subVolcano(contrFit2, lfc = lfc, padj = padj)
        }}

    if (length(cont.gen) == 2) {
        message(paste0("contrast is ",
                       paste(names(cont.gen), collapse = " vs ")))
        attr(cont.gen[[1]], "contrast") <- names(cont.gen)
        return(cont.gen[[1]])
    } else {
    return(cont.gen)
    }
    }
