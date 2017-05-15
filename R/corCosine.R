#' cosine correlation
#' @description Calculates the cosine correlation(s) between two vectors or
#' conformable matrices.
#' @export
#' @param x a numeric vector or matrix.
#' @param y a vector or matrix with compatible dimensions to x.
#' @details For two vectors, \eqn{X} and \eqn{Y}, the cosine similarity is
#' defined as \deqn{similarity=cos(\theta)=(\sum XY) \div (\sqrt(\sum X^2)\times \sqrt(\sum Y^2))}.
#' For centered data the cosine and Pearson's correlation coefficients are
#' equivalent.
#' @return a numeric matrix of cosine correlations.
#' @note \itemize{
#' \item{does not support missing values (`NA` nor `NaN`)}}
#' @seealso \code{\link{ntp}}
#' @references van Dongen S, Enright AJ. Metric distances derived from cosine
#' similarity and Pearson and Spearman correlations. arXiv:1208.3145 [cs, stat]
#' [Internet]. 2012 [cited 2016 Apr 22]; \url{http://arxiv.org/abs/1208.3145}
#' @examples
#' x <- stats::rnorm(1000, mean=0, sd=1)
#' y <- rbinom(1000, 1, 0.25)
#' corCosine(10*x, -x) # equals -1
#' corCosine(x, x)  # equals 1
#' replicate(10, corCosine(x, sample(x))) # expectation 0
#' hist(replicate(1000, corCosine(y, sample(x)))) # expectation 0
corCosine <- function(x, y) {
   x  <- as.matrix(x);y <- as.matrix(y)
    crossprod(x,y) /
            outer(sqrt(apply(x, 2, crossprod)), sqrt(apply(y, 2, crossprod)))
}
