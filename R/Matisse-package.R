#' Matisse: Multi-modal Analysis of Transcript Isoforms in Single-Cell
#' Sequencing Experiments
#'
#' @description
#' Matisse provides an integrated framework for isoform-resolved single-cell
#' RNA-seq analysis, built on top of \pkg{Seurat} and \pkg{Signac}.
#'
#' Key capabilities:
#' \itemize{
#'   \item \strong{MatisseObject} -- an S4 class that wraps a \code{Seurat}
#'     object and co-stores raw isoform counts, PSI matrices, and splice event
#'     annotations, keeping gene expression and isoform layers synchronised.
#'   \item \strong{PSI calculation} -- \code{\link{CalculatePSI}} computes
#'     per-cell Percent Spliced In values from raw junction or transcript
#'     counts. Works in both junction mode (STARsolo) and transcript mode
#'     (Bagpiper, FLAMES, LIQA).
#'   \item \strong{Isoform QC} -- per-cell metrics (\code{nCount_isoform},
#'     \code{nFeature_isoform}, \code{nPercent_isoform}) are written
#'     automatically at construction and by \code{\link{CalculatePSI}};
#'     \code{\link{FilterCells}} and \code{\link{FilterEvents}} enforce
#'     quality thresholds.
#'   \item \strong{Visualization} -- UMAP overlays, violin plots, PSI
#'     heatmaps, and sashimi junction-arc plots via a consistent
#'     ggplot2-based API.
#' }
#'
#' @section Package website:
#' Full documentation and vignettes are available at
#' \url{https://avisrilab.github.io/Matisse}.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom methods new is isVirtualClass validObject show setClass
#'   setGeneric setMethod setValidity slot initialize prototype as
#' @importFrom Matrix sparseMatrix nnzero rowSums colSums t Matrix
#' @importFrom rlang abort warn inform .data
#' @importFrom cli cli_alert_info cli_alert_success cli_alert_warning
#' @importFrom utils packageVersion
#' @importFrom stats median sd var hclust dist
## usethis namespace: end
NULL
