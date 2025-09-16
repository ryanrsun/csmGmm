#' symm_fit_cor_noAssumption.R
#'
#' Fit the correlated csmGmm for sets of correlated elements,
#' but we don't assume that the means in the composite alternative are greater in magnitude than those in the composite null.
#'
#' @param testStats J*K matrix of test statistics where J is the number of sets and K is number of elements in each set.
#' @param corMat K*K matrix that describes the correlation structure of each set.
#' @param initMuList List of 2^K elements where each element is a matrix with K rows and number of columns equal to the number of possible mean vectors for that binary group.
#' @param initPiList List of 2^K elements where each element is a vector with number of elements equal to the number of possible mean vectors for that binary group.
#' @param eps Scalar, stop the EM algorithm when L2 norm of difference in parameters is less than this value.
#' @param checkpoint Boolean, set to TRUE to print iterations of EM.
#'
#' @return A list with the elements:
#' \item{muInfo}{List with same dimensions as initMuList, holds the final mean parameters.}
#' \item{piInfo}{List with same dimensions as initPiList, holds the final probability parameters.}
#' \item{iter}{Number of iterations run in EM algorithm.}
#' \item{lfdrResults}{J*1 vector of all lfdr statistics.}
#' @importFrom mvtnorm rmvnorm
#' @importFrom dplyr %>% mutate arrange filter select slice
#' @import utils
#'
#' @export
#' @examples
#' set.seed(0)
#' corMat <- matrix(data=c(1, 0.3, 0.3, 1), nrow=2)
#' testStats <- rbind(mvtnorm::rmvnorm(n=200, mean=c(3, 0), sigma=corMat),
#' mvtnorm::rmvnorm(n=200, mean=c(0, 4), sigma=corMat),
#' mvtnorm::rmvnorm(n=100, mean=c(7, 7), sigma=corMat),
#' mvtnorm::rmvnorm(n=10^4 - 500, mean=c(0, 0), sigma=corMat))
#' initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3), nrow=2),
#' matrix(data=c(4, 0), nrow=2), matrix(data=c(5, 5), nrow=2))
#' initPiList <- list(c(0.9), c(0.04), c(0.04), c(0.02))
#' results <- symm_fit_cor_EM_noAssumption(testStats = testStats,
#' corMat = cor(testStats), initMuList = initMuList, initPiList = initPiList)
#'

symm_fit_cor_EM_noAssumption <- function(testStats, corMat, initMuList, initPiList, eps = 10^(-5), checkpoint=TRUE) {

  sigInv = solve(corMat)
  J <- nrow(testStats)
  # number of dimensions
  K <- ncol(testStats)
  B <- 2^K - 1
  # number of hl configurations - 1
  L <- 3^K - 1
  # make all configurations
  Hmat <- expand.grid(rep(list(-1:1), K))
  # attach the bl
  blVec <- rep(0, nrow(Hmat))
  slVec <- rep(0, nrow(Hmat))
  for (k_it in K:1) {
    blVec <- blVec + 2^(K - k_it) * abs(Hmat[, k_it])
    slVec <- slVec + abs(Hmat[, k_it])
  }
  # sort Hmat
  Hmat <- Hmat %>% dplyr::mutate(bl = blVec) %>%
    dplyr::mutate(sl = slVec) %>%
    dplyr::arrange(.data$bl, .data$Var1, .data$Var2) %>%
    dplyr::mutate(l = 0:(nrow(.) - 1))

  # initialize
  # the muInfo and piInfo are lists that hold the information in a more compact manner.
  # the allMu and allPi are matrices that repeat the information so the calculations can be performed faster.
  muInfo <- initMuList
  piInfo <- initPiList
  oldParams <- c(unlist(piInfo), unlist(muInfo))
  MbVec <- sapply(piInfo, FUN=length)

  # run until convergence
  diffParams <- 10
  iter <- 0
  while (diffParams > eps) {

    #####################################################
    # first, update allPi and allMu with the new parameter values

    # allPi holds the probabilities of each configuration (c 1), the bl of that
    # configuration (c 2), the sl of that configuration (c3), the m of that configuration (c 4),
    # and the l of that configuration (c 5),

    # number of rows in piMat is L if Mb = 1 for all b.
    # number of rows is \sum(l = 0 to L-1) {Mbl}
    allPi <- c(piInfo[[1]], 0, 0, 1, 0)

    # allMu holds the mean vectors for each configuration in each row, number of columns is K
    allMu <- rep(0, K)

    # loop through possible values of bl
    for (b_it in 1:B) {

      # Hmat rows with this bl
      tempH <- Hmat %>% dplyr::filter(.data$bl == b_it)

      # loop through possible m for this value of bl
      for (m_it in 1:MbVec[b_it + 1]) {
        allPi <- rbind(allPi, cbind(rep(piInfo[[b_it + 1]][m_it] / (2^tempH$sl[1]), nrow(tempH)), tempH$bl,
                                    tempH$sl, rep(m_it, nrow(tempH)), tempH$l))

        for (h_it in 1:nrow(tempH)) {
          allMu <- cbind(allMu, unlist(tempH %>% dplyr::select(-.data$bl, -.data$sl, -.data$l) %>%
                                         dplyr::slice(h_it)) * muInfo[[b_it + 1]][, m_it])
        }
      } # done looping through different m
    } # dont looping through bl
    colnames(allPi) <- c("Prob", "bl", "sl", "m", "l")

    ##############################################################################################
    # this is the E step where we calculate Pr(Z|c_l,m) for all c_l,m.
    # each l,m is one column.
    # only independence case for now
    conditionalMat <- sapply(X=data.frame(allMu), FUN=calc_dens_cor, Zmat = testStats, corMat = corMat) %>%
      sweep(., MARGIN=2, STATS=allPi[, 1], FUN="*")
    probZ <- apply(conditionalMat, 1, sum)
    AikMat <- conditionalMat %>% sweep(., MARGIN=1, STATS=probZ, FUN="/")
    Aik_alln <- apply(AikMat, 2, sum) / J

    ###############################################################################################
    # this is the M step for probabilities of hypothesis space
    for (b_it in 0:B) {
      for (m_it in 1:MbVec[b_it + 1]) {
        tempIdx <- which(allPi[, 2] == b_it & allPi[, 4] == m_it)
        piInfo[[b_it + 1]][m_it] <- sum(Aik_alln[tempIdx])
      }
    }

    # M step for the means
    # loop through values of bl
    # do the alternative last, enforce that it must be larger in magnitude than the nulls
    for (b_it in 1:B) {

      tempHmat <- Hmat %>% dplyr::filter(.data$bl == b_it)
      # loop through m
      for (m_it in 1:MbVec[b_it + 1]) {
        tempRightSum <- rep(0, nrow(allMu))
        tempLeftSum <- rep(0, nrow(allMu))

        # these are the classes that contribute to \mu_bl,m
        AikIdx <- which(allPi[, 2] == b_it & allPi[, 4] == m_it)
        for (idx_it in 1:length(AikIdx)) {
          tempAik <- AikIdx[idx_it]
          tempHvec <- tempHmat %>% dplyr::select(-.data$bl, -.data$sl, -.data$l) %>%
            dplyr::slice(idx_it) %>% unlist(.)
          LsigInvL <- diag(tempHvec) %*% sigInv %*% diag(tempHvec)

          tempLeftSum <- tempLeftSum + colSums(AikMat[, tempAik] * sweep(x = testStats %*% sigInv, MARGIN = 2,
                                                                         STATS = tempHvec, FUN="*"))
          tempRightSum <- tempRightSum + (J * Aik_alln[tempAik]) * LsigInvL
        } # done looping for one l, m

        # matrix inverse only for alternative
        if (b_it == B) {
          muInfo[[b_it + 1]][, m_it] <- solve(tempRightSum) %*% tempLeftSum
        } else {
          whichZero <- which(diag(tempRightSum) == 0)
          tempRightSum[whichZero, whichZero] <- 1
          muInfo[[b_it + 1]][, m_it] <- tempLeftSum / diag(tempRightSum)
        }

        # make sure mean constraint is satisfied
        #if (b_it == B) {
        #  maxMeans <- find_max_means(muInfo)
        #  whichSmaller <- which(muInfo[[b_it + 1]][, m_it] < maxMeans)
        #  if (length(whichSmaller) > 0) {
        #    muInfo[[b_it + 1]][whichSmaller, m_it] <- maxMeans[whichSmaller]
        #  }
        #} # done with mean constraint

      } # done looping through m

    } # done updating means

    ###############################################################################################
    # find difference
    allParams <- c(unlist(piInfo), unlist(muInfo))
    diffParams <- sum((allParams - oldParams)^2)

    # update
    oldParams <- allParams
    iter <- iter + 1
    if (checkpoint) {
      cat(iter, " - ", diffParams, "\n", allParams, "\n")
    }
  }

  # calculate local fdrs
  nullCols <- which(allPi[, 3] < K)
  probNull <- apply(conditionalMat[, nullCols], 1, sum)
  lfdrResults <- probNull / probZ

  return(list(piInfo = piInfo, muInfo = muInfo, iter = iter,
              lfdrResults = lfdrResults))
}
