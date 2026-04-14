#' Matisse: Multi-modal Analysis of Transcript Isoforms in Single-Cell
#' Sequencing Experiments
#'
#' @description
#' Matisse provides an integrated framework for isoform-resolved single-cell
#' RNA-seq analysis, built on top of \pkg{Seurat} and \pkg{Signac}.
#'
#' Key capabilities:
#' \itemize{
#'   \item \strong{MatisseObject} — an S4 class that wraps a \code{Seurat}
#'     object and co-stores junction counts, PSI matrices, and splice event
#'     annotations, keeping gene expression and isoform layers synchronised.
#'   \item \strong{PSI calculation} — \code{\link{CalculatePSI}} computes
#'     per-cell Percent Spliced In values from raw junction read counts and
#'     a user-supplied or auto-generated event annotation table.
#'   \item \strong{Isoform QC} — \code{\link{ComputeIsoformQC}} derives
#'     per-cell metrics (junction detection rate, event coverage, mean PSI);
#'     \code{\link{FilterCells}} and \code{\link{FilterEvents}} enforce
#'     quality thresholds.
#'   \item \strong{Visualization} — UMAP overlays, violin plots, PSI
#'     heatmaps, and junction coverage bar charts via a consistent
#'     ggplot2-based API.
#' }
#'
#' @section Package website:
#' Full documentation and vignettes are available at
#' \url{https://k3yavi.github.io/Matisse}.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom methods new is validObject show setClass setGeneric setMethod
#'   setValidity slot initialize prototype isVirtualClass
#' @importFrom Matrix sparseMatrix nnzero rowSums colSums t Matrix
#' @importFrom rlang abort warn inform check_installed
#' @importFrom cli cli_alert_info cli_alert_success cli_alert_warning
#'   cli_progress_bar cli_progress_update cli_progress_done
## usethis namespace: end
NULL
