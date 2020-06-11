#' nearest template prediction
#' @description Nearest Template Prediction (NTP) based on predefined class
#' templates.
#' @export
#' @param emat a numeric matrix with row features and sample columns.
#' \code{rownames(emat)} are matched against \code{templates$probe}.
#' @param templates a data frame with two columns; \emph{class} (coerced to
#' factor) and \emph{probe} (coerced to character).
#' @param nPerm an integer, number of permutations for \eqn{p}-value
#' estimation.
#' @param seed an integer, for \eqn{p}-value reproducibility. Setting seed
#' enforces serial processing.
#' @param distance a character, either c("cosine", "pearson", "spearman" or
#' "kendall").
#' @param nCores an integer specifying number of threads for
#' parallelization.
#' @param verbose logical, whether console messages are to be displayed.
#' @param doPlot logical, whether to produce prediction \code{\link{subHeatmap}}.
#' @details \code{ntp} implements the Nearest Template Prediction (NTP)
#' algorithm, largely as proposed by Yujin Hoshida (2010) (see below). For each
#' sample, distances to templates are calculated and class assigned based on
#' smallest distance. Distances are transformed from the sample-templates
#' correlations as follows:
#' \deqn{d.class = \sqrt(1/2 * (1-(cor(sample,templates))}
#' Template values are 1 for class features and 0 for non-class features (-1 if
#' there are only two classes). Prediction confidence is estimated based on
#' the distance of the null-distribution, estimated from permutation tests.
#' Thus the lowest possible estimate of the \eqn{p}-value is \eqn{1/nPerm}.
#' \itemize{
#'  \item{\code{emat}}{ should be a row-wise \emph{centered and scaled} matrix.
#'   For large, balanced datasets, this may be achieved by applying
#'   \code{\link{ematAdjust}} function.}
#'  \item{\code{templates}}{ is a data.frame defining class templates. A class
#'   template is a set of marker genes with higher expected expression in
#'   samples belonging to class compared to non-class samples. \code{templates}
#'   must contain at least two columns named \emph{probe} and \emph{class}.}
#'  \item{compared to Hoshida (2010), resulting \eqn{p}-value estimates are more
#'  conservative (by a factor equaling the number of classes) and
#'  the distances are a monotonic transformation of \eqn{1-cor} (see
#'  Details section above).}
#'  \item{Hoshida (2010) does not explicitly state whether input should be
#'  log2-transformed or not and examples includes both. Based on experience
#'  this choice affects results only at the margins, but for high-quality
#'  datasets, normalized, untransformed inputs may yield a small increase in
#'  accuracy.}
#'  }
#' For further details on the NTP algorithm, please refer to package vignette
#' and \href{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0015543}{Hoshida (2010)}.
#'
#' Parallel processing is implemented through \pkg{parallel}
#' \code{\link[parallel]{mclapply}} or \pkg{snow} \code{\link[snow]{parLapply}}
#' for nix and Windows systems, respectively.
#' @note \itemize{
#' \item{features with missing values are discarded.}
#' \item{setting seed disables parallel processing to ensure \eqn{p}-value
#' reproducibility.}
#' \item{for two random uncorrelated vectors \eqn{x,y} \eqn{N\sim(0,1)}
#' \eqn{E[d.xy]\approx0.71} when distance is cosine.}
#' \item{internally, correlations instead of distances are calculated.}
#' \item{accepts reuse of features (marker not specific for one class only)}}
#' @return a data frame with class predictions, template distances,
#' \eqn{p}-values and false discovery rate adjusted \eqn{p}-values
#' (\code{\link[stats]{p.adjust}(method = "fdr")}). Rownames equal emat
#' colnames.
#' @seealso \code{\link{corCosine}}, \code{\link[stats]{cor}}
#' @references Hoshida, Y. (2010). Nearest Template Prediction: A
#' Single-Sample-Based Flexible Class Prediction with Confidence Assessment.
#' PLoS ONE 5, e15543.
#' @examples
#' emat <- ematAdjust(crcTCGAsubset, normMethod = "quantile")
#' res <- ntp(emat, templates.CMS, doPlot=TRUE, nPerm=100)
#' head(res)
ntp <- function(emat, templates, nPerm = 1000, distance = "cosine",
                nCores = 1, seed = NULL, verbose = getOption("verbose"),
                doPlot = FALSE) {

    # checkInput ##############################################################

    # emat tests
    if (class(emat)[1] == "ExpressionSet")
        emat <- suppressPackageStartupMessages(Biobase::exprs(emat))
    if (is.data.frame(emat) | is.vector(emat)) emat <- as.matrix(emat)
    if (is.null(rownames(emat))) stop("missing emat rownames, check input")

    # templates tests
    if (is.null(templates$class) | is.null(templates$probe)) {
        stop("missing columns in templates, check input") }
    if (is.character(templates$class))
        templates$class <- factor(templates$class)
    if (is.factor(templates$probe) | is.integer(templates$probe)) {
        warning ("templates$probe coerced to character", call.=FALSE)
        templates$probe <-as.character(templates$probe)
    }

    if (is.character(distance) & isTRUE(verbose)) {
            message(paste0(distance, " correlation distance"))
    }

    # clean emat - for distCosine 0-imputation should be tested
    keepP <- stats::complete.cases(emat)
    if (sum(!keepP) > 0) {
        if (isTRUE(verbose)) message(paste0(sum(!keepP), "/", length(keepP),
                                            " features NA, discarded"))
        emat <- emat[keepP,,drop = FALSE]
    }

    # clean templates
    keepT <- templates$probe %in% rownames(emat)
    if (sum(!keepT) > 0) {
        if (isTRUE(verbose)) message(paste0(sum(!keepT), "/", length(keepT),
                " templates features not in emat, discarded"))
        templates <- templates[keepT,]
    }

    if (min(table(templates$class ))<2) {
        message ("<2 matched features/class")
        stop("check templates$probe is matchable against rownames(emat)",
                call. = FALSE)
    }

    if (min(table(templates$class ))<5) {
        warning("<5 matched features/class - unstable predictions",
                call.= FALSE)
    }

    # prepareInput ############################################################

    N <- ncol(emat)
    K <- nlevels(templates$class)
    S <- nrow(templates)
    P <- nrow(emat)

    class.names <- levels(templates$class)
    templates$class <- as.numeric(templates$class)

    # provide warning if emat seems non-normalized
    emat.mean <- round(mean(emat),2)
    if (abs(emat.mean) >1) {
        isnorm <- " <- check feature centering!"
        emat.sd <- round(stats::sd(emat),2)
        warning(paste0("emat mean=", emat.mean, "; sd=", emat.sd, isnorm),
                call.=FALSE)
    }

    # output classification overview
    feat.class <- paste(range(table(templates$class)),collapse = "-")
    if (isTRUE(verbose)) message(paste0(N, " samples; ",
                                        K, " classes; ",
                                        feat.class, " features/class"))

    # matching vector for emat and templates
    mm <- match(templates$probe, rownames(emat),nomatch = 0)

    if (!all(rownames(emat)[mm] == templates$probe)) {
        stop("error matching probes, check rownames(emat) and templates$probe")
    }

    # if features are reused across classes sample(..., replace=TRUE)
    pReplace <- length(templates$probe) > length(unique(templates$probe))

    # prepareTemplates ########################################################

    tmat <- matrix(rep(templates$class,K), ncol = K) # templates matrix
    for (k in seq_len(K)) tmat[,k] <- as.numeric(tmat[,k] == k)
    if (K == 2) tmat[tmat==0] <- -1

    # selectDistance ##########################################################

    if (!distance %in%
        c("cosine", "pearson", "spearman", "kendall"))
        stop("invalid distance method")

    if (distance == "cosine") {
        simFun <- function(x,y) corCosine(x,y)
    } else {
        simFun <- function(x, y) {
            stats::cor(x, y, method = distance)
        }
    }

    # ntpFunction #############################################################
    ntpFUN <- function(n) {

        # sample-templates correlations
        n.sim <- as.vector(simFun(emat[mm,n, drop = FALSE],tmat))

        # optimized for speed not readability
        # matrix(emat[,n][sample.int... makes permuted matrix
        # apply(simFun... calculates correlation and return max value

        n.sim.perm.max <- apply(simFun(
                matrix(emat[,n][sample.int(P, S*nPerm, replace=TRUE)],
                       ncol = nPerm), tmat), 1, max)


        n.ntp <- which.max(n.sim)
        # estimate p-value
        n.sim.ranks <- rank(-c(n.sim[n.ntp],(n.sim.perm.max)))
        n.pval <- n.sim.ranks[1]/length(n.sim.ranks)

        # return results
        return(c(
            n.ntp,                # prediction
            simToDist(n.sim),     # distance to all templates
            n.pval))              # p-value
    }

    # checkParalellization ####################################################
    # enforce serial processing if seed is set
    if (!is.null(seed)) {
        set.seed(seed)
        nCores <- 1
    }

    # check for paralellization dependencies
    existParallel <- packageExists("parallel")
    existSnow <- packageExists("snow")

    # paralellizedPrediction ##################################################
    # try parallelization if package is available and nCores not set to 1
    if ( (existParallel | existParallel) & nCores != 1) {

        funVal <- vector(mode = "numeric", length = 2+K)

    if ( .Platform$OS.type == "windows" ) {
            # win and *nix implementation using snow
            nSlaves <- ifelse(nCores == 0,
                              as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')),
                              nCores)
            if (nSlaves == 0) nSlaves <- parallel::detectCores()
            if (isTRUE(verbose)) message(paste0("parallel; ",
                                            nSlaves, " sockets pkg:snow",
                                            "; ", nPerm, " permutation(s)..."))

            # avoid un-neccessary dispatches
            nParts <- split(seq_len(N), cut(seq_len(N), nSlaves, labels=FALSE))

            cl <- snow::makeCluster(nSlaves, type = "SOCK")
            res <- snow::parLapply(cl, nParts, function(n)
                                            vapply(n, ntpFUN, funVal))
            snow::stopCluster(cl)
            res <- data.frame(t(do.call(cbind, res)))

            } else {
            # *nix implementation using parallel and mclapply
            # undocumented - user should specify nCores
            nCores <-  ifelse(nCores == 0,
                                parallel::detectCores(),
                                nCores)

            if (isTRUE(verbose)) message(paste0("parallel; ",
                                            nCores, " cores pkg:parallel",
                                            "; ", nPerm, " permutation(s)..."))
            options(mc.cores = nCores)

            # avoids un-neccessary dispatches but
            # memory and system.time scales with nPerm
            nParts <- split(seq_len(N), cut(seq_len(N), nCores, labels=FALSE))

            res <- parallel::mclapply(nParts, function(n)
                vapply(n, ntpFUN, funVal))
            res <- data.frame(t(do.call(cbind, res)))
            }
        } else {

        if (isTRUE(verbose)) message(paste0("serial processing; ",
                                            nPerm, " permutation(s)..."))

    # serializedPrediction ####################################################
        res <- lapply (seq_len(N), ntpFUN)
        res <- data.frame(do.call(rbind,res))
    }

    # prepareOutput ###########################################################
    colnames(res) <- c("prediction",
                       paste0("d.",class.names),
                      "p.value")
    res$prediction <- factor(class.names[res$prediction], levels = class.names)
    rownames(res) <- colnames(emat)
    res$p.value[res$p.value < 1/nPerm] <- 1/nPerm
    res$FDR <- stats::p.adjust(res$p.value, "fdr")

    # visualize ###############################################################
    if (isTRUE(doPlot)) {
        subHeatmap(emat = emat, res = res, templates = templates)
    }

    # returnOutput ############################################################
    if (isTRUE(verbose)) {
        message("predicted samples/class (FDR<0.05)")
        print(table(suppressMessages(subSetNA(res, FDR=0.05))$prediction,
                    useNA = "always"))
    }

    return(res)
}
