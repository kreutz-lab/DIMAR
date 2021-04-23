#' dimarPlotFacetHist
#'
#' @description Plots a faceted histogram displaying the percentage of missing values for
#' separate intensity intervals.
#' For this the median protein intensity values calculated across all samples are used.
#' # Addapted from:
# Mingyi Liu, Ashok Dongre, Proper imputation of missing values in proteomics datasets
# for differential expression analysis, Briefings in Bioinformatics, 2020;, bbaa112,
# https://doi.org/10.1093/bib/bbaa112
#' @return None
#' @param mtx Numeric matrix
#' @param nbins Amount of intensity intervals the range of intensity values should be split into.
#' @param savePlot Logical. If TRUE, the histogram is saved as pdf to working directory.
#' @param width Width of pdf to which plot is saved if savePlot = TRUE.
#' @param height Height of pdf to which plot is saved if savePlot = TRUE.
#' @export dimarPlotFacetHist
#' @examples
#' library(DIMAR)
#' library(openxlsx)
#' filename <- "Test2.xlsx"
#' filepath <- system.file("extdata", filename, package = "DIMAR")
#' df <- read.xlsx(filepath, sheet="Sheet1", startRow = 2)
#' row.names(df) <- paste(c(1:nrow(df)), df$`Protein.(Uniprot.ID)`, sep = "_")
#' df <- df[, grepl("^AD\\d|^C\\d", names(df))]
#' mtx <- as.matrix(df)
#' dimarPlotFacetHist(mtx)

dimarPlotFacetHist <- function(mtx, nbins = 20, savePlot = FALSE, width = 13, height = 5) {
  df <- as.data.frame(mtx)
  hist.df <- data.frame(Na.perc = apply(df, 1, function(row)sum(is.na(row))/length(row)*100),
                        medianIntensities = apply(df, 1, stats::median, na.rm = TRUE))
  # Bins of equal length
  hist.df$medianIntensities_i <- ggplot2::cut_interval(hist.df$medianIntensities, nbins)

  histplot <- ggplot2::ggplot(hist.df, ggplot2::aes(Na.perc, group = medianIntensities_i)) +
    ggplot2::geom_histogram(binwidth = 1) +
    ggplot2::labs(x = paste0("Missing % over ", ncol(df), " samples in ", nbins, " intensity intervals"), y = "Count") +
    ggplot2::facet_wrap(~ medianIntensities_i, nrow = 1) +
    ggplot2::theme_gray() +
    ggplot2::theme(strip.text.x = ggplot2::element_text(size = 8),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
  if (savePlot) {
    ggplot2::ggsave(histplot, filename = "facethistplot.pdf", width = width, height = height)

  } else {
    histplot
  }
}
