#' dimarPlotHeatmap
#'
#' @description Plots a heatmap of the protein intensity of a dataset. The samples are displayed
#' in the columns and the proteins in the rows. The proteins are sorted according to an increasing
#' number of missing values and a decreasing mean across all samples.
#' @return None
#' @param mtx Numeric matrix
#' @param savePlot Logical. If TRUE, the histogram is saved as pdf to working directory.
#' @param width Width of pdf to which plot is saved if savePlot = TRUE.
#' @param height Height of pdf to which plot is saved if savePlot = TRUE.
#' @export dimarPlotHeatmap
#' @examples
#' library(DIMAR)
#' library(openxlsx)
#' filename <- "Test2.xlsx"
#' filepath <- system.file("extdata", filename, package = "DIMAR")
#' df <- read.xlsx(filepath, sheet="Sheet1", startRow = 2)
#' row.names(df) <- paste(c(1:nrow(df)), df$`Protein.(Uniprot.ID)`, sep = "_")
#' df <- df[, grepl("^AD\\d|^C\\d", names(df))]
#' mtx <- as.matrix(df)
#' dimarPlotHeatmap(mtx)

dimarPlotHeatmap <- function(mtx, savePlot = FALSE, width = 6, height = 6) {

  mtx.heatmap <- mtx[order(-rowSums(is.na(mtx)), rowMeans(mtx, na.rm = TRUE)), ]

  if (savePlot) {
    grDevices::pdf("heatmap.pdf", width = width, height = height)
    heatmap(mtx.heatmap, Rowv = NA, Colv = NA, na.rm = T,
            col = RColorBrewer::brewer.pal(n = 9, name = "Blues"), scale = "none",
            labRow = FALSE, margins = c(6, 0.5))
    grDevices::dev.off()
  } else {
    heatmap(mtx.heatmap, Rowv = NA, Colv = NA, na.rm = T,
            col = RColorBrewer::brewer.pal(n = 9, name = "Blues"), scale = "none",
            labRow = FALSE, margins = c(6, 0.5))
  }
}
