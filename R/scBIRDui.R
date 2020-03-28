#' scBIRDui
#'
#' Launch the BIRD user interface in local machine
#'
#' This function will automatically launch the SCRAT user interface in a web browser.
#' The user interface can also be accessed by xxx. Neither R nor any packages are required in this online version.
#' However, it is highly recommended that the user interface be launched locally for faster running speed.
#'
#' @export
#' @import shiny shinyBS GenomicAlignments ggplot2 reshape2 pheatmap scatterD3 gplots DT mclust tsne dbscan
#' @author Yueyi Li, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' \dontrun{
#'    scBIRDui()
#' }

scBIRDui <- function() {
  shiny::runApp(system.file("shiny",package="scBIRD"))
}
