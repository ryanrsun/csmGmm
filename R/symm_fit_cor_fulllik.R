#' symm_fit_cor_fulllik.R
#'
#' Full likelihood, block correlation, blocks of size 2
#'
#' @param testStats J*K matrix of test statistics where J is the number of sets and K is number of elements in each set.
#' @param corMat K*K matrix that describes the correlation structure of each 2 by 2 block.
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
#'
#' @importFrom dplyr %>% filter select slice
#' @importFrom magrittr %>% set_colnames
#' @import utils
#'
#' @export
#' @examples
#' set.seed(0)
#' testStats <- cbind(rnorm(10^4), rnorm(10^4))
#' testStats[1:100, 1] <- rnorm(100, mean=3)
#' testStats[101:200, 1] <- rnorm(100, mean=5)
#' testStats[201:300, 2] <- rnorm(100, mean=4)
#' testStats[301:400, 1:2] <- rnorm(200, mean=7)
#' initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3), nrow=2, ncol=1),
#' matrix(data=c(3, 0), nrow=2, ncol=1), matrix(data=c(6, 6), nrow=2, ncol=1))
#' initPiList <- list(c(0.9), c(0.04),c(0.04), c(0.02))
#' results <- symm_fit_cor_EM_fulllik(testStats = testStats, corMat=diag(c(1,1)),
#' initMuList = initMuList, initPiList = initPiList)
#'

symm_fit_cor_EM_fulllik <- function(testStats, corMat, initMuList, initPiList, eps = 10^(-5), checkpoint=TRUE) {

  # need to split the test statistics for block correlation
  datK1 <- matrix(data=testStats[, 1], ncol=2, byrow=TRUE)
  datK2 <- matrix(data=testStats[, 2], ncol=2, byrow=TRUE)

  # rho
  rho <- corMat[1, 2]

  # number of composite hypotheses
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
  Hmat <- Hmat %>% mutate(bl = blVec) %>%
    mutate(sl = slVec) %>%
    arrange(.data$bl, .data$Var1, .data$Var2) %>%
    mutate(l = 0:(nrow(.) - 1))

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

    # For full likelihood block correlation approach, we need 9*9=81 possible mu vectors
    # It's really only 25 (since can be 0, -1, 1, -3, 3 for each dimension), but this makes code easier
    # R expands sets of 1:9, 1:9 by going 1:9 of the first column, all with 1 of the second column,
    # then 1:9 of the first column, all with 2 of the second column, etc.
    # This is the opposite of the way I want, but since my two columns are equal, I can just swap
    # the order of the two columns of the end result.
    # Then transpose to be like allMu, so I can do the sapply correctly
    muK1 = t(expand.grid(allMu[1,], allMu[1,])[, c(2, 1)])
    muK2 = t(expand.grid(allMu[2,], allMu[2,])[, c(2, 1)])
    muIdxMat <- expand.grid(1:9, 1:9)[, c(2, 1)] %>%
      magrittr::set_colnames(c("K1", "K2"))
    twoDPi <- apply(expand.grid(allPi[, 1], allPi[, 1]), 1, prod)

    ##############################################################################################
    # this is the E step where we calculate Pr(Z|c_l,m) for all c_l,m.
    # each l,m is one column.
    # only independence case for now
    conditionalMatK1log <- sapply(X = data.frame(muK1), FUN=calc_dens_cor, Zmat=datK1, corMat=corMat, log=TRUE) %>%
      sweep(., MARGIN=2, STATS=log(twoDPi), FUN="+")
    conditionalMatK2log <- sapply(X = data.frame(muK2), FUN=calc_dens_cor, Zmat=datK2, corMat=corMat, log=TRUE) %>%
      sweep(., MARGIN=2, STATS=twoDPi, FUN="+")
    # these are numerator of the B_a_l,l'
    BnumeratorLog <- conditionalMatK1log + conditionalMatK2log
    Bnumerator <- exp(BnumeratorLog)
    BDenom <- apply(Bnumerator, 1, sum)
    # the B_a_l,l' (what used to be the AikMat)
    BllMatlog <- sweep(BnumeratorLog, MARGIN=1, STATS=log(BDenom), FUN="-")
    # we're just a little inaccurate here, not sure what to do with tiny probabilities...
    BllMat <- exp(BllMatlog)

    # the marginals (one half of the Bll, the Aik type terms)
    sumEachNine <- function(x) {
      c(sum(x[1:9]), sum(x[10:18]), sum(x[19:27]), sum(x[28:36]), sum(x[37:45]),
        sum(x[46:54]), sum(x[55:63]), sum(x[64:72]), sum(x[73:81]))
    }
    sumEveryNine <- function(x) {
      c(sum(x[c(1,10,19,28,37,46,55,64,73)]), sum(x[c(2,11,20,29,38,47,56,65,74)]),
        sum(x[c(3,12,21,30,39,48,57,66,75)]), sum(x[c(4,13,22,31,40,49,58,67,76)]),
        sum(x[c(5,14,23,32,41,50,59,68,77)]), sum(x[c(6,15,24,33,42,51,60,69,78)]),
        sum(x[c(7,16,25,34,43,52,61,70,79)]), sum(x[c(8,17,26,35,44,53,62,71,80)]),
        sum(x[c(9,18,27,36,45,54,63,72,81)]))
    }
    BllMat_alln <- apply(BllMat, 2, sum) / (J / 2)
    # the odd j
    AikMat1 <- t(apply(BllMat, 1, sumEachNine))
    Aik1_alln <- apply(AikMat1, 2, sum) / J
    # the even j
    AikMat2 <- t(apply(BllMat, 1, sumEveryNine))
    Aik2_alln <- apply(AikMat2, 2, sum) / J

    ###############################################################################################
    # this is the M step for probabilities of hypothesis space
    for (b_it in 0:B) {
      for (m_it in 1:MbVec[b_it + 1]) {
        tempIdx <- which(allPi[, 2] == b_it & allPi[, 4] == m_it)
        piInfo[[b_it + 1]][m_it] <- sum(Aik1_alln[tempIdx] + Aik2_alln[tempIdx])
      }
    }

    ###############################################################################################
    # M step for means
    # equations for first dimension mean mu_2,1
    eq21term21 <- 0
    eq21term31 <- 0
    eq21term0 <- 0

    # which l=1-9 holds the means I care about
    AikIdx21 <- which(allPi[, 2] == 2 & allPi[, 4] == 1)
    # forward is it's in the first element
    # reverse it mu_2,1 is in the second element
    colAikIdx21 <- which(muIdxMat$K1 %in% AikIdx21 | muIdxMat$K2 %in% AikIdx21 )

    # loop through
    testDF <- data.frame(a1 = rep(0, length(colAikIdx21)), a2=0, a3=0, a4=0, a5=0, a6=0)
    for (idx_it in 1:length(colAikIdx21)) {
      # columns of BllMat
      tempBllCol <- colAikIdx21[idx_it]

      # figure out which of the three cases I'm in
      tempHmatRow <- unlist(muIdxMat[tempBllCol, ])
      tempBl1 <- Hmat$bl[tempHmatRow[1]]
      tempBl2 <- Hmat$bl[tempHmatRow[2]]
      # one mu_2,1 and one 0
      if (tempBl1 == 2 & tempBl2 %in% c(0, 1)) {
        tempL <- sign(Hmat$Var1[tempHmatRow[1]])
        eq21term21 <- eq21term21 - sum(BllMat[, tempBllCol]) / (1 - rho^2)
        eq21term0 <- eq21term0 + sum(BllMat[, tempBllCol] * tempL * (datK1[, 1] - rho * datK1[, 2])) / (1 - rho^2)
      }
      # one 0 and one mu_2,1
      if (tempBl1 %in% c(0, 1) & tempBl2 == 2) {
        # remember we're in the first dimension here so it should still be var1
        tempL <- sign(Hmat$Var1[tempHmatRow[2]])
        eq21term21 <- eq21term21 - sum(BllMat[, tempBllCol]) / (1 - rho^2)
        eq21term0 <- eq21term0 + sum(BllMat[, tempBllCol] * tempL * (datK1[, 2] - rho * datK1[, 1])) / (1 - rho^2)
      }
      # one mu_2,1 and one mu_3,1
      if ((tempBl1 == 2 & tempBl2 == 3) | (tempBl1 == 3 & tempBl2 == 2)) {
        tempL1 <- sign(Hmat$Var1[tempHmatRow[1]])
        tempL2 <- sign(Hmat$Var1[tempHmatRow[2]])
        tempLcross <- tempL1 * tempL2
        eq21term21 = eq21term21 - sum(BllMat[, tempBllCol]) / (1 - rho^2)
        eq21term31 = eq21term31 + sum(BllMat[, tempBllCol] * rho * tempLcross) / (1 - rho^2)
        if (tempBl1 == 1 & tempBl2 == 3) {
          eq21term0 <- eq21term0 + sum(BllMat[, tempBllCol] * tempL1 * (datK1[, 1] - rho * datK1[, 2])) / (1 - rho^2)
        } else {
          eq21term0 <- eq21term0 + sum(BllMat[, tempBllCol] * tempL2 * (datK1[, 2] - rho * datK1[, 1])) / (1 - rho^2)
        }
      }
      # two mu_2,1
      if (tempBl1 == 2 & tempBl2 == 2) {
        tempL1 <- sign(Hmat$Var1[tempHmatRow[1]])
        tempL2 <- sign(Hmat$Var1[tempHmatRow[2]])
        tempLcross <- tempL1 * tempL2
        eq21term21 <- eq21term21 + sum(BllMat[, tempBllCol] * (-2 + 2 * tempLcross * rho)) / (1 - rho^2)
        eq21term0 <- eq21term0 + sum(BllMat[, tempBllCol] * (tempL1 * (datK1[, 1] - rho * datK1[, 2]) + tempL2 * (datK1[, 2] - rho * datK1[, 1]))) / (1 - rho^2)
      }
      #cat(c(eq21term21, eq21term0))
      testDF[idx_it, ] <- c(eq21term21, eq21term0, BllMat_alln[tempBllCol], muK1[, tempBllCol], twoDPi[tempBllCol])
    } # done looping through all idx for mu_2,1

    #-------------------------------------------------------#
    # equations for first dimension mean mu_3,1

    # which l=1-9 holds the means I care about
    AikIdx31 <- which(allPi[, 2] == 3 & allPi[, 4] == 1)
    # forward is it's in the first element
    # reverse it mu_2,1 is in the second element
    colAikIdx31 <- which(muIdxMat$K1 %in% AikIdx31 | muIdxMat$K2 %in% AikIdx31)

    eq31term21 <- 0
    eq31term31 <- 0
    eq31term0 <- 0
    # loop through
    for (idx_it in 1:length(colAikIdx31)) {
      # columns of BllMat
      tempBllCol <- colAikIdx31[idx_it]

      # figure out which of the three cases I'm in
      tempHmatRow <- unlist(muIdxMat[tempBllCol, ])
      tempBl1 <- Hmat$bl[tempHmatRow[1]]
      tempBl2 <- Hmat$bl[tempHmatRow[2]]
      # one mu_3,1 and one 0
      if (tempBl1 == 3 & tempBl2 %in% c(0, 1)) {
        tempL <- sign(Hmat$Var1[tempHmatRow[1]])
        eq31term31 <- eq31term31 - sum(BllMat[, tempBllCol]) / (1 - rho^2)
        eq31term0 <- eq31term0 + sum(BllMat[, tempBllCol] * tempL * (datK1[, 1] - rho * datK1[, 2])) / (1 - rho^2)
      }
      # one 0 and one mu_3,1
      if (tempBl1 %in% c(0, 1) & tempBl2 == 3) {
        tempL <- sign(Hmat$Var1[tempHmatRow[2]])
        eq31term31 <- eq31term31 - sum(BllMat[, tempBllCol]) / (1 - rho^2)
        eq31term0 <- eq31term0 + sum(BllMat[, tempBllCol] * tempL * (datK1[, 2] - rho * datK1[, 1])) / (1 - rho^2)
      }
      # one mu_3,1 and one mu_2,1
      if ((tempBl1 == 3 & tempBl2 == 2) | (tempBl1 == 2 & tempBl2 == 3)) {
        tempL1 <- sign(Hmat$Var1[tempHmatRow[1]])
        tempL2 <- sign(Hmat$Var1[tempHmatRow[2]])
        tempLcross <- tempL1 * tempL2
        eq31term31 = eq31term31 - sum(BllMat[, tempBllCol]) / (1 - rho^2)
        eq31term21 = eq31term21 + sum(BllMat[, tempBllCol] * rho * tempLcross) / (1 - rho^2)
        if (tempBl1 == 3 & tempBl2 == 1) {
          eq31term0 <- eq31term0 + sum(BllMat[, tempBllCol] * tempL1 * (datK1[, 1] - rho * datK1[, 2])) / (1 - rho^2)
        } else {
          eq31term0 <- eq31term0 + sum(BllMat[, tempBllCol] * tempL2 * (datK1[, 2] - rho * datK1[, 1])) / (1 - rho^2)
        }
      }
      # two mu_3,1
      if (tempBl1 == 3 & tempBl2 == 3) {
        tempL1 <- sign(Hmat$Var1[tempHmatRow[1]])
        tempL2 <- sign(Hmat$Var1[tempHmatRow[2]])
        tempLcross <- tempL1 * tempL2
        eq31term31 <- eq31term31 + sum(BllMat[, tempBllCol] * (-2 + 2 * tempLcross * rho)) / (1 - rho^2)
        eq31term0 <- eq31term0 + sum(BllMat[, tempBllCol] * (tempL1 * (datK1[, 1] - rho * datK1[, 2]) + tempL2 * (datK1[, 2] - rho * datK1[, 1]))) / (1 - rho^2)
      }
    } # done looping through all idx for mu_3,1

    # now solve the equations
    Amat <- matrix(data=c(eq21term21, eq21term31, eq31term21, eq31term31), nrow=2, byrow=T)
    cVec <- matrix(data=c(eq21term0, eq31term0), nrow=2)
    muCol1 <- solve(Amat) %*% -cVec

    ###############################################################################################
    # now do the second set of means
    # equations for second dimension mean mu_1,2
    eq12term12 <- 0
    eq12term32 <- 0
    eq12term0 <- 0

    # which l=1-9 holds the means I care about
    AikIdx12 <- which(allPi[, 2] == 1 & allPi[, 4] == 1)
    # forward is it's in the first element
    # reverse it is in the second element
    colAikIdx12 <- which(muIdxMat$K1 %in% AikIdx12 | muIdxMat$K2 %in% AikIdx12)

    # loop through
    for (idx_it in 1:length(colAikIdx12)) {
      # columns of BllMat
      tempBllCol <- colAikIdx12[idx_it]

      # figure out which of the three cases I'm in
      tempHmatRow <- unlist(muIdxMat[tempBllCol, ])
      tempBl1 <- Hmat$bl[tempHmatRow[1]]
      tempBl2 <- Hmat$bl[tempHmatRow[2]]
      # one mu_1,2 and one 0
      if (tempBl1 == 1 & tempBl2 %in% c(0, 2)) {
        tempL <- sign(Hmat$Var2[tempHmatRow[1]])
        eq12term12 <- eq12term12 - sum(BllMat[, tempBllCol]) / (1 - rho^2)
        eq12term0 <- eq12term0 + sum(BllMat[, tempBllCol] * tempL * (datK2[, 1] - rho * datK2[, 2])) / (1 - rho^2)
      }
      # one 0 and one mu_1,2
      if (tempBl1 %in% c(0, 2) & tempBl2 == 1) {
        tempL <- sign(Hmat$Var2[tempHmatRow[2]])
        eq12term12 <- eq12term12 - sum(BllMat[, tempBllCol]) / (1 - rho^2)
        eq12term0 <- eq12term0 + sum(BllMat[, tempBllCol] * tempL * (datK2[, 2] - rho * datK2[, 1])) / (1 - rho^2)
      }
      # one mu_1,2 and one mu_3,2
      if ((tempBl1 == 1 & tempBl2 == 3) | (tempBl1 == 3 & tempBl2 == 1)) {
        tempL1 <- sign(Hmat$Var2[tempHmatRow[1]])
        tempL2 <- sign(Hmat$Var2[tempHmatRow[2]])
        tempLcross <- tempL1 * tempL2
        eq12term12 = eq12term12 - sum(BllMat[, tempBllCol]) / (1 - rho^2)
        eq12term32 = eq12term32 + sum(BllMat[, tempBllCol] * rho * tempLcross) / (1 - rho^2)
        if (tempBl1 == 1 & tempBl2 == 3) {
          eq12term0 <- eq12term0 + sum(BllMat[, tempBllCol] * tempL1 * (datK2[, 1] - rho * datK2[, 2])) / (1 - rho^2)
        } else {
          eq12term0 <- eq12term0 + sum(BllMat[, tempBllCol] * tempL2 * (datK2[, 2] - rho * datK2[, 1])) / (1 - rho^2)
        }
      }
      # two mu_1,2
      if (tempBl1 == 1 & tempBl2 == 1) {
        tempL1 <- sign(Hmat$Var2[tempHmatRow[1]])
        tempL2 <- sign(Hmat$Var2[tempHmatRow[2]])
        tempLcross <- tempL1 * tempL2
        eq12term12 <- eq12term12 + sum(BllMat[, tempBllCol] * (-2 + 2 * tempLcross * rho)) / (1 - rho^2)
        eq12term0 <- eq12term0 + sum(BllMat[, tempBllCol] * (tempL1 * (datK2[, 1] - rho * datK2[, 2]) + tempL2 * (datK2[, 2] - rho * datK2[, 1]))) / (1 - rho^2)
      }
    } # done looping through all idx for mu_1,2

    #-------------------------------------------------------#
    # equations for second dimension mean mu_3,2
    eq32term12 <- 0
    eq32term32 <- 0
    eq32term0 <- 0

    # which l=1-9 holds the means I care about
    AikIdx32 <- which(allPi[, 2] == 3 & allPi[, 4] == 1)
    # forward is it's in the first element
    # reverse it mu_2,1 is in the second element
    colAikIdx32 <- which(muIdxMat$K1 %in% AikIdx32 | muIdxMat$K2 %in% AikIdx32)

    # loop through
    for (idx_it in 1:length(colAikIdx32)) {
      # columns of BllMat
      tempBllCol <- colAikIdx32[idx_it]

      # figure out which of the three cases I'm in
      tempHmatRow <- unlist(muIdxMat[tempBllCol, ])
      tempBl1 <- Hmat$bl[tempHmatRow[1]]
      tempBl2 <- Hmat$bl[tempHmatRow[2]]
      # one mu_3,2 and one 0
      if (tempBl1 == 3 & tempBl2 %in% c(0, 1)) {
        tempL <- sign(Hmat$Var2[tempHmatRow[1]])
        eq32term32 <- eq32term32 - sum(BllMat[, tempBllCol]) / (1 - rho^2)
        eq32term0 <- eq32term0 + sum(BllMat[, tempBllCol] * tempL * (datK2[, 1] - rho * datK2[, 2])) / (1 - rho^2)
      }
      # one 0 and one mu_3,2
      if (tempBl1 %in% c(0, 1) & tempBl2 == 3) {
        tempL <- sign(Hmat$Var2[tempHmatRow[2]])
        eq32term32 <- eq32term32 - sum(BllMat[, tempBllCol]) / (1 - rho^2)
        eq32term0 <- eq32term0 + sum(BllMat[, tempBllCol] * tempL * (datK2[, 2] - rho * datK2[, 1])) / (1 - rho^2)
      }
      # one mu_3,2 and one mu_1,2
      if ((tempBl1 == 3 & tempBl2 == 1) | (tempBl1 == 1 & tempBl2 == 3)) {
        tempL1 <- sign(Hmat$Var2[tempHmatRow[1]])
        tempL2 <- sign(Hmat$Var2[tempHmatRow[2]])
        tempLcross <- tempL1 * tempL2
        eq32term32 = eq32term32 - sum(BllMat[, tempBllCol]) / (1 - rho^2)
        eq32term12 = eq32term12 + sum(BllMat[, tempBllCol] * rho * tempLcross) / (1 - rho^2)
        if (tempBl1 == 3 & tempBl2 == 1) {
          eq32term0 <- eq32term0 + sum(BllMat[, tempBllCol] * tempL1 * (datK2[, 1] - rho * datK2[, 2])) / (1 - rho^2)
        } else {
          eq32term0 <- eq32term0 + sum(BllMat[, tempBllCol] * tempL2 * (datK2[, 2] - rho * datK2[, 1])) / (1 - rho^2)
        }
      }
      # two mu_3,2
      if (tempBl1 == 3 & tempBl2 == 3) {
        tempL1 <- sign(Hmat$Var2[tempHmatRow[1]])
        tempL2 <- sign(Hmat$Var2[tempHmatRow[2]])
        tempLcross <- tempL1 * tempL2
        eq32term32 <- eq32term32 + sum(BllMat[, tempBllCol] * (-2 + 2 * tempLcross * rho)) / (1 - rho^2)
        eq32term0 <- eq32term0 + sum(BllMat[, tempBllCol] * (tempL1 * (datK2[, 1] - rho * datK2[, 2]) + tempL2 * (datK2[, 2] - rho * datK2[, 1]))) / (1 - rho^2)
      }
    } # done looping through all idx for mu_3,2

    # now solve the equations
    Amat2 <- matrix(data=c(eq12term12, eq12term32, eq32term12, eq32term32), nrow=2, byrow=T)
    cVec2 <- matrix(data=c(eq12term0, eq32term0), nrow=2)
    muCol2 <- solve(Amat2) %*% -cVec2

    # update muInfo
    # muCol1 is for mu_2,1 and mu_3,1
    # muCol2 is for mu_1,2 and mu_3,2
    muInfo[[2]][2,1] <- muCol2[1]
    muInfo[[2]][1,1] <- 0
    muInfo[[3]][1,1] <- muCol1[1]
    muInfo[[3]][2,1] <- 0
    muInfo[[4]][1,1] <- muCol1[2]
    muInfo[[4]][2,1] <- muCol2[2]

    # mean constraint
    maxMeans <- find_max_means(muInfo)
    whichSmaller <- which(muInfo[[4]][, 1] < maxMeans)
    if (length(whichSmaller) > 0) {
      muInfo[[4]][whichSmaller, 1] <- maxMeans[whichSmaller]
    } # done with mean constraint

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
  } # done with EM

  # calculate local fdrs
  nullCols <- which(allPi[, 3] < K)
  probNull1 <- apply(AikMat1[, nullCols], 1, sum)
  probNull2 <- apply(AikMat2[, nullCols], 1, sum)
  rearranged <- cbind(probNull1, probNull2)
  lfdrResults <- as.numeric(t(rearranged))
  #lfdrResults <- probNull / probZ

  return(list(piInfo = piInfo, muInfo = muInfo, iter = iter,
              lfdrResults = lfdrResults))
}
