#' create_plots.R
#'
#' Creates lfdr-value Manhttan plots and z-score plots to visualize composite null hypothesis testing results
#'
#' @param plotData data.frame with columns labeled (exactly) Chr, BP, cum_avg_lfdr, Z_Dataset1, Z_Dataset2
#'
#' @return A list with the elements:
#' \item{manPlot}{A ggplot object of the Manhattan plot of the lfdr-values at which each SNP is significant}
#' \item{zPlot}{A ggplot object of a scatterplot comparing the z-scores for the first two datasets}
#'
#' @import dplyr ggplot2
#' @importFrom rlang .data
#'
#' @export
#' @examples
#' set.seed(0)
#' plotData <- data.frame(Chr = rep(1:22, each=100), BP=runif(n=2200, min=1, max=10^8),
#' cum_avg_lfdr=runif(n=2200, min=0, max=22),
#' Z_Dataset1=rnorm(n=2200), Z_Dataset2=rnorm(n=2200))
#' create_plots(plotData = plotData)
#'
create_plots <- function(
    plotData
) {

  # -------------------------------
  # Format data
  # -------------------------------
  plotData <- plotData %>%
    mutate(
      Chr = as.numeric(.data$Chr),
      BP = as.numeric(.data$BP),
      cum_avg_lfdr = as.numeric(.data$cum_avg_lfdr)
    ) %>%
    filter(
      !is.na(.data$Chr),
      !is.na(.data$BP),
      !is.na(.data$cum_avg_lfdr)
    )

  # -------------------------------
  # Chromosome lengths
  # -------------------------------
  chrCounts <- plotData %>%
    dplyr::group_by(.data$Chr) %>%
    dplyr::summarise(
      chrLength = max(.data$BP, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$Chr) %>%
    dplyr::pull(.data$chrLength)

  # -------------------------------
  # Create manhattan plot
  # -------------------------------

  # arrange data by chromosome
  plotData <- plotData %>% arrange(.data$Chr)
  uniqueChrs <- sort(unique(plotData$Chr))

  chrOffsets <- cumsum(chrCounts)
  truePos <- rep(NA, nrow(plotData))
  counter <- 1
  for (chr_it in 1:length(uniqueChrs)) {
    tempChr <- uniqueChrs[chr_it]
    tempDat <- plotData %>% filter(.data$Chr == tempChr)
    # offset = sum of all previous chromosomes' max positions
    offsetVal <- ifelse(tempChr == 1, 0, chrOffsets[tempChr - 1])
    truePos[counter:(counter + nrow(tempDat) - 1)] <- offsetVal + tempDat$BP
    counter <- counter + nrow(tempDat)
  }

  # x-axis ticks at the end of each chromosome
  xBreaks <- chrOffsets[1:22]
  xBreaksLabs <- ifelse((1:22) %% 2 == 0, "", as.character(1:22))

  plotData <- plotData %>% mutate(truePos = truePos)

  manPlot <- ggplot(plotData, aes(x=.data$truePos, y=-log10(.data$cum_avg_lfdr))) +
    geom_point(color = "blue", shape = 16) +
    xlab("Chromosome") + ylab("-log10(lfdr)") +
    scale_x_continuous(name="Chr", breaks=xBreaks, labels=xBreaksLabs) +
    coord_cartesian(ylim = c(0,max(-log10(plotData$cum_avg_lfdr)) + 1)) +
    theme_bw()

  # -------------------------------
  # Create Z-score plot
  # -------------------------------
  zPlot <- ggplot(plotData, aes(x = .data$Z_Dataset1, y = .data$Z_Dataset2)
  ) +
    geom_point(color = "red", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    xlab(paste0("Z-score 1")) +
    ylab(paste0("Z-score 2")) +
    theme_bw(base_size = 16)

  return(list(
    manPlot = manPlot,
    zPlot = zPlot))
}

