#' calc_dens_ind.R
#'
#' Calculate the density of bivariate normal under fitted acsGmm.
#'
#' @param x 2*1 vector of means.
#' @param Zmat J*2 matrix of test statistics.
#'
#' @return A J*1 vector of densities for each row of Zmat.
#'
#' @export
#' @examples
#' x <- c(0, 0)
#' Zmat <- cbind(rnorm(10^5), rnorm(10^5))
#' calc_dens_ind_2d(x, Zmat)
#'
calc_dens_ind_2d <- function(x, Zmat) {
  # right now, only for 2 or 3 dimensions
  dnorm(Zmat[, 1], mean=x[1], sd=1) * dnorm(Zmat[, 2], mean=x[2], sd=1)
}


#' Calculate the density of trivariate normal under fitted acsGmm.
#'
#' @param x 3*1 vector of means.
#' @param Zmat J*3 matrix of test statistics.
#'
#' @return A J*1 vector of densities for each row of Zmat.
#'
#' @export
#' @examples
#' x <- c(0, 0)
#' Zmat <- cbind(rnorm(10^5), rnorm(10^5), rnorm(10^5))
#' calc_dens_ind_3d(x, Zmat)
#'
calc_dens_ind_3d <- function(x, Zmat) {
  # right now, only for 2 or 3 dimensions
  dnorm(Zmat[, 1], mean=x[1], sd=1) * dnorm(Zmat[, 2], mean=x[2], sd=1) * dnorm(Zmat[, 3], mean=x[3], sd=1)
}



#' Calculate the density of multivariate normal under fitted acsGmm.
#'
#' @param x K*1 vector of means.
#' @param Zmat J*K matrix of test statistics.
#'
#' @return A J*1 vector of densities for each row of Zmat.
#'
#' @export
#' @examples
#' x <- c(0, 0)
#' Zmat <- cbind(rnorm(10^5), rnorm(10^5), rnorm(10^5), rnorm(10^5))
#' calc_dens_ind_multiple(x, Zmat)
#'
calc_dens_ind_multiple <- function(x, Zmat) {
  K <- ncol(zMat)
  tempProd <- rep(1, nrow(zMat))
  for (k_it in 1:K) {
    tempProd <- tempProd * dnorm(Zmat[, k_it], mean=x[k_it], sd=1)
  }
}
