voomTransform <- function(emat, normMethod, design=NULL ,simplify=FALSE, ...) {

        validMethods <- c(subData$methods.voom, subData$methods.edgeR)
        if (!normMethod %in% validMethods) stop("invalid normalization method")

        # courtsey of Roman Lustrik http://stackoverflow.com/questions/3476782/
        is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)
            all(abs(as.vector(x) - round(as.vector(x))) < tol)

        if (!is.wholenumber(emat)) warning ("emat is not integers")

        if (!packageExists("limma"))
            stop (paste0("voom requires limma package available from:\n",
                         "http://www.bioconductor.org/packages/release/bioc/html/limma.html"))


        if (normMethod %in% subData$methods.edgeR) {
            if (!packageExists("edgeR"))
                stop (paste0("normalization method requires edgeR package available from:\n",
                             "http://www.bioconductor.org/packages/release/bioc/html/edgeR.html"))
            dge <- edgeR::DGEList(counts = emat)
            dge <- edgeR::calcNormFactors(dge, method = normMethod, ...)
            emat <- limma::voom(dge, design, plot = FALSE)
        }

        if (normMethod %in% subData$methods.voom) {
            emat <- limma::voom(emat, design, plot = FALSE,
                                normalize.method = normMethod)
        }

    if (simplify==TRUE) {
        P <- ncol(emat$E)
        rowN <- rownames(emat$E)
        emat <- emat$E
        rownames(emat) <- rowN

    }
    return(emat)
    }
