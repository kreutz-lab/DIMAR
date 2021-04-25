#' dimarPlotSampleDensity
#'
#' @description Plots the protein intensity distributions of each sample.
#' @return None
#' @param mtx Numeric matrix
#' @param savePlot Logical. If TRUE, the histogram is saved as pdf to working directory.
#' @param width Width of pdf to which plot is saved if savePlot = TRUE.
#' @param height Height of pdf to which plot is saved if savePlot = TRUE.
#' @export dimarPlotSampleDensity
#' @examples
#' library(DIMAR)
#' library(openxlsx)
#' filename <- "Test2.xlsx"
#' filepath <- system.file("extdata", filename, package = "DIMAR")
#' df <- read.xlsx(filepath, sheet="Sheet1", startRow = 2)
#' row.names(df) <- paste(c(1:nrow(df)), df$`Protein.(Uniprot.ID)`, sep = "_")
#' df <- df[, grepl("^AD\\d|^C\\d", names(df))]
#' mtx <- as.matrix(df)
#' dimarPlotSampleDensity(mtx)

dimarPlotSampleDensity <- function(mtx, savePlot = FALSE, width = 5, height = 13) {
  df <- as.data.frame(mtx)
  df.long <- utils::stack(df)
  colnames(df.long) <- c("Intensity", "Sample")

  ridgeplot <- ggplot2::ggplot(df.long, ggplot2::aes(x = Intensity, y = forcats::fct_rev(Sample))) +
    ggridges::geom_density_ridges(alpha = 0) +
    ggplot2::theme_minimal()

  if (savePlot) {
    ggplot2::ggsave(ridgeplot, filename = "densityplot.pdf", width = width, height = height)
  } else {
    ridgeplot
  }
}
