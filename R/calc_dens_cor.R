#' calc_dens_cor.R
#'
#' For J*2 matrix of J sets, calculate the density of bivariate normal under fitted c-csmGmm.
#'
#' @param x 2*1 vector of means.
#' @param Zmat J*2 matrix of test statistics.
#' @param corMat 2*2 matrix describing correlation structure of test statistics.
#' @param log return logarithm of density
#'@importFrom mvtnorm dmvnorm
#'
#' @return A J*1 vector of densities for each row of Zmat.
#'
#' @export
#' @examples
#' x <- c(0, 0)
#' Zmat <- cbind(rnorm(10^5), rnorm(10^5))
#' calc_dens_cor(x, Zmat, corMat = cor(Zmat))
#'
calc_dens_cor <- function(x, Zmat, corMat, log=FALSE) {
  # right now, only for 2 or 3 dimensions
  # will do each row of Zmat separately
  mvtnorm::dmvnorm(Zmat, mean=x, sigma = corMat, log=log)
}
