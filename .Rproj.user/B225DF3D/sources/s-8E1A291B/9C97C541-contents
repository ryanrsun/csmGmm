#' plot_findings.R
#'
#' Plot findings (Z1 vs. Z2, colored based on significant or not), from a single output of fitting).
#'
#' @param testStatistics J*K matrix of test statistics where J is the number of sets and K is number of elements in each set.
#' @param lfdrVec J*1 vector of all lfdr statistics.
#' @param cutoff Cutoff for significance.
#' @param numPts Number of points to actually plot, if plotting all may take a long time to render. Randomly sample points.
#' @param plotPDF Boolean, if TRUE then color significant according to lfdr value, not average (cumulative mean) Lfdr of region.
#' @param xlimits 2*1 vector of limits on x-axis.
#' @param ylimits 2*1 vector of limits on y-axis.
#' @param plotPts Vector of points to plot, in case want to match them up across different figures.
#'
#' @return A ggplot
#'
#' @export
#' @examples
#' set.seed(0)
#' testStatistics <- cbind(rnorm(10^5), rnorm(10^5))
#' lfdrVec <- runif(n=10^5)
#' plot_findings <- function(testStatistics = testStatistics, lfdrVec = lfdrVec, cutoff = 0.1,
#' numPts = 20000, plotPDF = TRUE, xlimits=c(-10, 10), ylimits=c(-10, 10), plotPts=NULL)
#'
plot_findings <- function(testStatistics, lfdrVec, cutoff, numPts = 20000, plotPDF = FALSE,
                          xlimits=c(-10, 10), ylimits=c(-10, 10), plotPts=NULL) {

  # arrange data
  plotDat <- testStatistics %>% as.data.frame(.) %>%
    set_colnames(paste0("Z", 1:ncol(testStatistics))) %>%
    mutate(lfdr = lfdrVec) %>%
    mutate(idx = 1:nrow(.)) %>%
    arrange(lfdr) %>%
    mutate(LfdrAvg = cummean(lfdr))

  if (!is.null(plotPts)) {
    plotDat <- plotDat %>% filter(idx %in% plotPts)
  } else {
    plotDat <- plotDat %>% slice_sample(n=numPts)
  }

  # significant by running average or straight localfdr?
  if (plotPDF) {
    plotDat <- plotDat %>% mutate(Sig = ifelse(lfdr < cutoff, 1, 0))
  } else {
    plotDat <- plotDat %>% mutate(Sig = ifelse(LfdrAvg < cutoff, 1, 0))
  }

  # plot
  returnPlot <- ggplot(data=plotDat, aes(x = Z1, y=Z2, color=as.factor(Sig))) +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    theme_cowplot() +
    xlim(xlimits) + ylim(ylimits) +
    scale_color_manual(name = "Significant", values = c("#F8766D", "#00BFC4"), labels = c("FALSE", "TRUE"))

  return(returnPlot)
}
