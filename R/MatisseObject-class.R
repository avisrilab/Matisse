#' The MatisseObject S4 class
#'
#' The central data structure for Matisse. It wraps a \code{\link[Seurat]{Seurat}}
#' object and augments it with isoform-resolved layers: raw junction counts, per-
#' event PSI matrices, and splice event annotations. All isoform matrices use the
#' same cell barcodes as the embedded Seurat object, keeping the two layers
#' synchronized.
#'
#' @slot seurat A \code{Seurat} object carrying gene-level expression,
#'   dimensionality reductions, and cell metadata.
#' @slot psi A sparse matrix (dgCMatrix, cells x events) of PSI values in
#'   \eqn{[0, 1]}. \code{NA} entries indicate insufficient read coverage.
#' @slot inclusion_counts A sparse matrix (dgCMatrix, cells x events) of
#'   aggregated inclusion-junction read counts per event.
#' @slot exclusion_counts A sparse matrix (dgCMatrix, cells x events) of
#'   aggregated exclusion-junction read counts per event.
#' @slot junction_counts A sparse matrix (dgCMatrix, cells x junctions) of raw
#'   per-junction read counts.
#' @slot event_data A \code{data.frame} with one row per splice event. Required
#'   columns: \code{event_id}, \code{gene_id}, \code{chr}, \code{strand},
#'   \code{event_type}, \code{inclusion_junctions}, \code{exclusion_junctions}.
#' @slot junction_data A \code{data.frame} with one row per junction. Required
#'   columns: \code{junction_id}, \code{chr}, \code{start}, \code{end},
#'   \code{strand}, \code{gene_id}.
#' @slot isoform_metadata A \code{data.frame} of per-cell isoform QC metrics.
#'   Rownames correspond to cell barcodes.
#' @slot version Character string recording the Matisse version used to create
#'   the object.
#' @slot misc Named list for user-defined extra data.
#'
#' @name MatisseObject-class
#' @rdname MatisseObject-class
#' @exportClass MatisseObject
setClass(
  "MatisseObject",
  slots = c(
    seurat            = "ANY",        # Seurat object
    psi               = "ANY",        # dgCMatrix | NULL
    inclusion_counts  = "ANY",        # dgCMatrix | NULL
    exclusion_counts  = "ANY",        # dgCMatrix | NULL
    junction_counts   = "ANY",        # dgCMatrix | NULL
    event_data        = "data.frame",
    junction_data     = "data.frame",
    isoform_metadata  = "data.frame",
    version           = "character",
    misc              = "list"
  ),
  prototype = list(
    seurat            = NULL,
    psi               = NULL,
    inclusion_counts  = NULL,
    exclusion_counts  = NULL,
    junction_counts   = NULL,
    event_data        = data.frame(),
    junction_data     = data.frame(),
    isoform_metadata  = data.frame(),
    version           = as.character(utils::packageVersion("Matisse")),
    misc              = list()
  )
)

# ---------------------------------------------------------------------------
# Validity
# ---------------------------------------------------------------------------

setValidity("MatisseObject", function(object) {
  errors <- character()

  # seurat slot must be a Seurat object if non-NULL
  if (!is.null(object@seurat)) {
    if (!inherits(object@seurat, "Seurat")) {
      errors <- c(errors,
        "'seurat' slot must be a Seurat object or NULL.")
    }
  }

  cells <- .get_cells(object)

  # If a PSI matrix exists, its row names must match cell barcodes
  if (!is.null(object@psi)) {
    if (!is.null(cells) && !identical(rownames(object@psi), cells)) {
      errors <- c(errors,
        "Row names of 'psi' must match cell barcodes in the Seurat object.")
    }
  }

  # inclusion_counts and exclusion_counts must have the same dimensions as psi
  if (!is.null(object@psi) && !is.null(object@inclusion_counts)) {
    if (!identical(dim(object@psi), dim(object@inclusion_counts))) {
      errors <- c(errors,
        "'inclusion_counts' must have the same dimensions as 'psi'.")
    }
    if (!identical(colnames(object@psi), colnames(object@inclusion_counts))) {
      errors <- c(errors,
        "Column names of 'inclusion_counts' must match those of 'psi'.")
    }
  }

  if (!is.null(object@psi) && !is.null(object@exclusion_counts)) {
    if (!identical(dim(object@psi), dim(object@exclusion_counts))) {
      errors <- c(errors,
        "'exclusion_counts' must have the same dimensions as 'psi'.")
    }
  }

  # event_data must have required columns if non-empty
  if (nrow(object@event_data) > 0) {
    required_cols <- c("event_id", "gene_id", "chr", "strand",
                       "event_type", "inclusion_junctions",
                       "exclusion_junctions")
    missing_cols <- setdiff(required_cols, colnames(object@event_data))
    if (length(missing_cols) > 0) {
      errors <- c(errors, paste0(
        "event_data is missing required columns: ",
        paste(missing_cols, collapse = ", ")))
    }
  }

  # junction_data must have required columns if non-empty
  if (nrow(object@junction_data) > 0) {
    required_cols <- c("junction_id", "chr", "start", "end",
                       "strand", "gene_id")
    missing_cols <- setdiff(required_cols, colnames(object@junction_data))
    if (length(missing_cols) > 0) {
      errors <- c(errors, paste0(
        "junction_data is missing required columns: ",
        paste(missing_cols, collapse = ", ")))
    }
  }

  if (length(errors) == 0) TRUE else errors
})

# ---------------------------------------------------------------------------
# Internal helpers used inside validity and elsewhere
# ---------------------------------------------------------------------------

# Extract cell barcodes from the embedded Seurat object (or NULL)
.get_cells <- function(object) {
  if (is.null(object@seurat)) return(NULL)
  colnames(object@seurat)
}

# Number of cells
.n_cells <- function(object) {
  cells <- .get_cells(object)
  if (is.null(cells)) 0L else length(cells)
}

# Number of splice events
.n_events <- function(object) {
  if (is.null(object@psi)) return(0L)
  ncol(object@psi)
}
