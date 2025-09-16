#' check_incongruous.R
#'
#' Check the number of sets of test statistics that have a higher (less significant) lfdr value
#' than other sets with test statistics of uniformly smaller magnitudes.
#'
#' @param zMatrix J*K vector of all test statistics.
#' @param lfdrVec J*1 vector of lfdr values corresponding to each set of test statistics.
#'
#' @return A vector with all the indices of all sets that have a higher lfdr value those a set
#' with smaller test statistic magnitudes.
#' @importFrom dplyr mutate arrange desc select
#' @importFrom rlang .data
#'
#' @export
#' @examples
#' zMatrix <- cbind(rnorm(10^4), rnorm(10^4))
#' lfdrVec <- runif(10^4)
#' check_incongruous(zMatrix = zMatrix, lfdrVec = lfdrVec)
#'
check_incongruous <- function(zMatrix, lfdrVec) {

  # remove lfdr = 1
  lessThanOne <- which(lfdrVec < 0.99)
  if (length(lessThanOne) <= 1) {return(c())}
  zMatrix <- zMatrix[lessThanOne, ]
  lfdrVec <- lfdrVec[lessThanOne]

  # do it in K^2 quadrants
  K <- ncol(zMatrix)
  quadrants <- expand.grid(rep(list(c(-1, 1)), K))

  badIdx <- c()
  for (quad_it in 1:nrow(quadrants)) {
    # separate into quadrants
    idxVec <- 1:nrow(zMatrix)
    tempStats <- zMatrix
    tempLfdr <- lfdrVec
    for (k_it in 1:K) {
      if (class(tempStats)[1] == "numeric") {break}
      if (quadrants[quad_it, k_it] == -1) {
        toKeep <- which(tempStats[, k_it] < 0 )
        idxVec <- idxVec[toKeep]
        tempLfdr <- tempLfdr[toKeep]
      } else {
        toKeep <- which(tempStats[, k_it] > 0 )
        idxVec <- idxVec[toKeep]
        tempLfdr <- tempLfdr[toKeep]
      }
      tempStats <- tempStats[toKeep, ]
    } # end finding quadrant
    if (length(idxVec) <= 1) {next}
    # take absolute value
    tempStats <- abs(tempStats)

    # order by lfdr
    tempDat <- tempStats %>% as.data.frame(.data) %>%
      dplyr::mutate(lfdr = tempLfdr) %>%
      dplyr::mutate(idx = idxVec)

    # check for incongruous
    if (K == 2) {
      colnames(tempDat)[1:2] <- c("Z1", "Z2")
      tempDat <- tempDat %>%
        dplyr::arrange(tempLfdr, dplyr::desc(.data$Z1), dplyr::desc(.data$Z2))
      incongruousVec <- sapply(1:nrow(tempDat),FUN = find_2d,  allTestStats = as.matrix(tempDat %>% dplyr::select(.data$Z1, .data$Z2)))
    } else if (K == 3) {
      colnames(tempDat)[1:3] <- c("Z1", "Z2", "Z3")
      tempDat <- tempDat %>%
        dplyr::arrange(tempLfdr, dplyr::desc(.data$Z1), dplyr::desc(.data$Z2), dplyr::desc(.data$Z3))
      incongruousVec <- sapply(1:nrow(tempDat), FUN = find_3d, allTestStats = as.matrix(tempDat %>% dplyr::select(.data$Z1, .data$Z2, .data$Z3)))
    } else {
      stop("only support for 2-3 dimensions right now")
    }

    # get the bad indices
    badIdx <- c(badIdx, tempDat$idx[which(incongruousVec > 0)])
  }

  return(badIdx)
}


#' Tells if row x if allTestStats is an incongruous result (has a higher lfdr than a set of
#' test statistics with lower magnitudes). For K=2 case.
#'
#' @param x Scalar, which row of allTestStats to check.
#' @param allTestStats J*K vector of all test statistics.
#'
#' @return A scalar denoting the number of sets with lower lfdr and test statistics of lower magnitude. 0 means congruous result.
#'
#' @export
#' @examples
#' zMatrix <- cbind(rnorm(10^4), rnorm(10^4))
#' find_2d(x = 5, allTestStats = zMatrix)
#'
find_2d <- function(x, allTestStats) {
  length(which(allTestStats[1:x, 1] < allTestStats[x, 1] & allTestStats[1:x, 2] < allTestStats[x, 2]))
}

#' Tells if row x if allTestStats is an incongruous result (has a higher lfdr than a set of
#' test statistics with lower magnitudes). For K=3 case.
#'
#' @param x Scalar, which row of allTestStats to check.
#' @param allTestStats J*K vector of all test statistics.
#'
#' @return A scalar denoting the number of sets with lower lfdr and test statistics of lower magnitude. 0 means congruous result.
#'
#' @export
#' @examples
#' zMatrix <- cbind(rnorm(10^4), rnorm(10^4),  rnorm(10^4))
#' find_3d(x = 5, allTestStats = zMatrix)
#'
find_3d <- function(x, allTestStats) {
  length(which(allTestStats[1:x, 1] < allTestStats[x, 1] & allTestStats[1:x, 2] < allTestStats[x, 2] &
                 allTestStats[1:x, 3] < allTestStats[x, 3]))
}

