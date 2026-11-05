#' diagnostic_plot.R
#'
#' Collapse a fitted csmGmm into the one-dimensional marginal mixture for each dimension.
#'
#' Under the independence model, the marginal distribution of column k of testStats is a
#' univariate Gaussian mixture with unit variance. Every configuration where dimension k is null
#' contributes mass at mean 0, and every configuration where dimension k is non-null contributes
#' half its mass at +mu and half at -mu.
#'
#' @param fitResults The list returned by symm_fit_ind_EM() (only muInfo and piInfo are used).
#'
#' @return A data.frame with columns dim (dimension index), comp ("null", "negAlt", or "posAlt"),
#' mean, and weight. Weights sum to 1 within each dimension.
#'
#' @importFrom dplyr %>% group_by summarise
#' @importFrom rlang .data
#' @export
#' @examples
#' initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3, 0, 5), nrow=2, ncol=2),
#' matrix(data=c(3, 0, 5, 0), nrow=2, ncol=2), matrix(data=c(7, 7), nrow=2, ncol=1))
#' initPiList <- list(c(0.9), c(0.02, 0.02),c(0.02, 0.02), c(0.02))
#' get_marginal_mixture(list(muInfo = initMuList, piInfo = initPiList))
#'
get_marginal_mixture <- function(fitResults) {

  muInfo <- fitResults$muInfo
  piInfo <- fitResults$piInfo
  K <- nrow(muInfo[[1]])
  B <- length(piInfo) - 1

  compList <- list()
  for (b_it in 0:B) {
    # dimensions that are non-null for this bl, column 1 is the most significant bit (matches Hmat)
    altDims <- bitwAnd(as.integer(b_it), as.integer(2^(K - 1:K))) > 0

    for (m_it in seq_along(piInfo[[b_it + 1]])) {
      tempPi <- piInfo[[b_it + 1]][m_it]
      for (k_it in 1:K) {
        if (altDims[k_it]) {
          tempMu <- muInfo[[b_it + 1]][k_it, m_it]
          compList[[length(compList) + 1]] <- data.frame(dim = k_it, comp = c("negAlt", "posAlt"),
                                                         mean = c(-tempMu, tempMu), weight = tempPi / 2)
        } else {
          compList[[length(compList) + 1]] <- data.frame(dim = k_it, comp = "null", mean = 0, weight = tempPi)
        }
      }
    }
  }

  do.call(rbind, compList) %>%
    dplyr::group_by(.data$dim, .data$comp, .data$mean) %>%
    dplyr::summarise(weight = sum(.data$weight), .groups = "drop") %>%
    as.data.frame()
}


#' Histogram of the observed test statistics with the fitted marginal mixture density overlaid.
#'
#' For each dimension, plots a histogram of the observed test statistics along with the fitted
#' marginal mixture density (scaled to counts), the null N(0,1) component, and the alternative
#' component (all non-null components combined). Pink bars show observed counts minus the counts expected under the fitted null
#' component, similar to the diagnostics in repfdr::ztobins().
#'
#' @param testStats J*K matrix of test statistics that was passed to symm_fit_ind_EM().
#' @param fitResults The list returned by symm_fit_ind_EM().
#' @param nBins Number of histogram bins.
#' @param zLim Optional length-2 vector of x-axis limits; defaults to the range of each column.
#' Statistics outside zLim are dropped from the histogram but still count toward the total J.
#' @param yLim Optional length-2 vector of y-axis limits. The alternative components are usually tiny
#' relative to the null, so setting e.g. yLim = c(0, 50) zooms in on the tails without dropping any bars.
#'
#' @return A list of K ggplot objects, one for each dimension.
#'
#' @import ggplot2
#' @importFrom dplyr %>% filter
#' @importFrom graphics hist
#' @importFrom stats dnorm pnorm
#' @importFrom rlang .data
#' @export
#' @examples
#' set.seed(0)
#' testStats <- cbind(rnorm(10^4), rnorm(10^4))
#' testStats[1:200, 1] <- rnorm(200, mean=3)
#' testStats[201:400, 1] <- rnorm(200, mean=5)
#' testStats[401:600, 2] <- rnorm(200, mean=3)
#' testStats[601:800, 2] <- rnorm(200, mean=5)
#' testStats[801:1000, 1:2] <- rnorm(400, mean=7)
#' initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3, 0, 5), nrow=2, ncol=2),
#' matrix(data=c(3, 0, 5, 0), nrow=2, ncol=2), matrix(data=c(7, 7), nrow=2, ncol=1))
#' initPiList <- list(c(0.9), c(0.02, 0.02),c(0.02, 0.02), c(0.02))
#' results <- symm_fit_ind_EM(testStats = testStats, sameDirAlt = TRUE,
#' initMuList = initMuList, initPiList = initPiList)
#' histPlots <- plot_fit_hist(testStats = testStats, fitResults = results)
#' histPlots[[1]]
#'
plot_fit_hist <- function(testStats, fitResults, nBins = 100, zLim = NULL, yLim = NULL) {

  J <- nrow(testStats)
  K <- ncol(testStats)
  allComps <- get_marginal_mixture(fitResults)

  compLabels <- c(total = "Fitted mixture", null = "Null N(0,1)", alt = "Alternative")
  compColors <- c("Fitted mixture" = "#1b9e3e", "Null N(0,1)" = "#2c5fbf", "Alternative" = "#d62728")
  compTypes <- c("Fitted mixture" = "solid", "Null N(0,1)" = "dotted", "Alternative" = "dashed")

  plotList <- list()
  for (k_it in 1:K) {
    tempZ <- testStats[, k_it]
    tempLim <- if (is.null(zLim)) range(tempZ) else zLim
    tempComps <- allComps %>% dplyr::filter(.data$dim == k_it)

    # observed counts
    breaks <- seq(tempLim[1], tempLim[2], length.out = nBins + 1)
    binWidth <- breaks[2] - breaks[1]
    counts <- graphics::hist(tempZ[tempZ >= tempLim[1] & tempZ <= tempLim[2]], breaks = breaks, plot = FALSE)$counts
    mids <- (breaks[-1] + breaks[-length(breaks)]) / 2

    # counts expected under the null component, using exact bin probabilities
    nullWeight <- sum(tempComps$weight[tempComps$comp == "null"])
    expNull <- J * nullWeight * diff(stats::pnorm(breaks))
    binDat <- rbind(data.frame(mid = mids, count = counts, bar = "Observed"),
                    data.frame(mid = mids, count = pmax(counts - expNull, 0), bar = "Observed - expected null"))

    # fitted densities on the count scale
    zGrid <- seq(tempLim[1], tempLim[2], length.out = 500)
    densMat <- sapply(1:nrow(tempComps), function(x) {
      tempComps$weight[x] * stats::dnorm(zGrid, mean = tempComps$mean[x])
    })
    densMat <- matrix(densMat, nrow = length(zGrid))
    curveDat <- data.frame(z = zGrid, count = J * binWidth * rowSums(densMat), comp = compLabels[["total"]])
    # one alternative density
    for (comp_it in c("null", "alt")) {
      whichComp <- if (comp_it == "null") which(tempComps$comp == "null") else which(tempComps$comp != "null")
      if (length(whichComp) == 0) {next}
      curveDat <- rbind(curveDat, data.frame(z = zGrid, comp = compLabels[[comp_it]],
                                             count = J * binWidth * rowSums(densMat[, whichComp, drop = FALSE])))
    }
    curveDat$comp <- factor(curveDat$comp, levels = unname(compLabels))

    plotList[[k_it]] <- ggplot() +
      geom_col(data = binDat, aes(x = .data$mid, y = .data$count, fill = .data$bar),
               width = binWidth, position = "identity", color = "grey35", linewidth = 0.2) +
      geom_line(data = curveDat, aes(x = .data$z, y = .data$count, color = .data$comp, linetype = .data$comp),
                linewidth = 0.9) +
      scale_fill_manual(values = c("Observed" = "#b3b6bc", "Observed - expected null" = "#f4a6c0"), name = NULL) +
      scale_color_manual(values = compColors, name = NULL, drop = TRUE) +
      scale_linetype_manual(values = compTypes, name = NULL, drop = TRUE) +
      coord_cartesian(ylim = yLim) +
      labs(x = "Test statistic", y = "Count",
           title = paste0("csmGmm fitted marginal density, dimension ", k_it),
           subtitle = paste0("Fitted marginal null proportion = ", signif(nullWeight, 4))) +
      theme_bw() +
      theme(legend.position = "inside", legend.position.inside = c(0.99, 0.99),
            legend.justification = c(1, 1), legend.background = element_rect(fill = "white", color = "grey70"),
            legend.box.just = "right", legend.spacing.y = unit(0, "pt"))
  }

  plotList
}
