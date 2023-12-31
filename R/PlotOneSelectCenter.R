#' Plot one selected ROI
#'
#' @param sce A single cell experiment object.
#' @param ploti A number of indicate which ROI to plot.
#' @param enhanced A logical variable for using enhanced data or not.
#'
#' @importFrom gridExtra grid.arrange
#'
#' @return A figure object for the selected ROI.
#' @export
#' @examples
#' data("example_sce")
#' example_sce <- mySpatialPreprocess(example_sce, platform="Visium")
#' PlotOneSelectedCenter(example_sce, ploti = 1)
#'
PlotOneSelectedCenter <- function(sce,
                                  ploti,
                                  enhanced = FALSE) {
  
    stopifnot(exprs = {
      is.numeric(ploti)
    })
  
    if(is.null(S4Vectors::metadata(sce)$selectCenters$selectEnt)) {
        p1 <- FindRegionalCells(sce,
                               centerID = S4Vectors::metadata(sce)$selectCenters$selectCenters[ploti],
                               radius = S4Vectors::metadata(sce)$selectCenters$selectRadius[ploti],
                               enhanced = enhanced,
                               avern = 5,
                               doPlot = FALSE,
                               returnPlot = TRUE)
        p2 <- p1 + ggplot2::ggtitle(paste0("Point ", ploti, " Radius", S4Vectors::metadata(sce)$selectCenters$selectRadius[ploti]))
    } else {
        p1 <- FindRegionalCells(sce,
                               centerID = S4Vectors::metadata(sce)$selectCenters$selectID[ploti],
                               radius = S4Vectors::metadata(sce)$selectCenters$selectRadius[ploti],
                               enhanced = enhanced,
                               avern = 5,
                               doPlot = FALSE,
                               returnPlot = TRUE)
        p2 <- p1 + ggplot2::ggtitle(paste0("Point ", ploti, " Radius", S4Vectors::metadata(sce)$selectCenters$selectRadius[ploti],
                                          " Entropy", signif(S4Vectors::metadata(sce)$selectCenters$selectEnt[ploti],4)))
    }

    return(p2)
}
