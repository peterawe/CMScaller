# set global color palette used for plot functions
.onLoad <- function (libname, pkgname) {
    options(subClassCol = subData$classCol)
}

# good practice to enable silencing
.onAttach <- function (libname, pkgname) {
    packageStartupMessage("CMScaller v0.9.2; GENCODE v26/GRCh38.p10 (Brainarray v22)")
}

# internal utility functions for consistent read and write
writeTab <- function(mat, file, ...) {
    utils::write.table(mat, file, quote = FALSE, sep = "\t",
            row.names = FALSE, ...)
    }

readTab <- function(file, ...) {
    utils::read.table(file, header = TRUE, sep = "\t", ...)
    }
# convert from/to distance/correlation. Squared term removes sign
distToSim  <- function(x) 1-2*x^2
simToDist <- function(x) sqrt(1/2*(1-(x)))

IPR <- function(x, probs = c(.1, .9))
    abs(diff(stats::quantile(x, probs, names = FALSE, na.rm = TRUE)))


# checks whether package is installed
packageExists <- function (pkg) {
    stopifnot(is.character(pkg))
    pkg %in% utils::installed.packages()[,"Package"]
}

accuracy <- function(x,y) cohensKappa(x,y, adjusted=FALSE)


# identify coordinates for margin lines
# courtsey of jbaums http://stackoverflow.com/questions/30765866/
line2user <- function(line, side) {
    lh <- graphics::par('cin')[2] * graphics::par('cex') * graphics::par('lheight')
    x_off <- diff(graphics::grconvertX(c(0, lh), 'inches', 'npc'))
    y_off <- diff(graphics::grconvertY(c(0, lh), 'inches', 'npc'))
    switch(side,
           `1` = graphics::grconvertY(-line * y_off, 'npc', 'user'),
           `2` = graphics::grconvertX(-line * x_off, 'npc', 'user'),
           `3` = graphics::grconvertY(1 + line * y_off, 'npc', 'user'),
           `4` = graphics::grconvertX(1 + line * x_off, 'npc', 'user'),
           stop("Side must be 1, 2, 3, or 4", call.=FALSE))
}
