#' The MatisseObject S4 class
#'
#' The central data structure for Matisse. It wraps a \code{\link[Seurat]{Seurat}}
#' object and augments it with isoform-resolved layers. Two fixed assays live
#' inside the embedded Seurat object:
#' \describe{
#'   \item{\code{"transcript"}}{A standard \code{Assay5} holding raw
#'     transcript-level counts (transcripts × cells). Created by
#'     \code{\link{CreateMatisseObjectFromTranscripts}}.}
#'   \item{\code{"psi"}}{A \code{ChromatinAssay} (Signac) holding PSI values
#'     in the \code{"data"} layer, inclusion counts in \code{"counts"}, and
#'     exclusion counts in \code{"exclusion"} (all features × cells). Created
#'     by \code{\link{CalculatePSI}} or
#'     \code{\link{CreateMatisseObjectFromTranscripts}}.}
#' }
#' Raw junction counts and splice-event annotations are kept as slots on the
#' MatisseObject itself.
#'
#' @slot seurat A \code{Seurat} object carrying gene-level expression,
#'   dimensionality reductions, cell metadata, and the \code{"transcript"} /
#'   \code{"psi"} assays.
#' @slot junction_counts A sparse matrix (dgCMatrix, cells × junctions) of raw
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
    seurat           = "ANY",        # Seurat object (contains "transcript" Assay5 and "psi" ChromatinAssay)
    junction_counts  = "ANY",        # dgCMatrix | NULL (cells x junctions)
    event_data       = "data.frame",
    junction_data    = "data.frame",
    isoform_metadata = "data.frame",
    version          = "character",
    misc             = "list"
  ),
  prototype = list(
    seurat           = NULL,
    junction_counts  = NULL,
    event_data       = data.frame(),
    junction_data    = data.frame(),
    isoform_metadata = data.frame(),
    version          = as.character(utils::packageVersion("Matisse")),
    misc             = list()
  )
)

# ---------------------------------------------------------------------------
# Validity
# ---------------------------------------------------------------------------

setValidity("MatisseObject", function(object) {
  errors <- character()

  if (!is.null(object@seurat)) {
    if (!inherits(object@seurat, "Seurat")) {
      errors <- c(errors, "'seurat' slot must be a Seurat object or NULL.")
    }
  }

  cells <- .get_cells(object)

  # If junction_counts exists, its rows must match cell barcodes
  if (!is.null(object@junction_counts)) {
    if (!is.null(cells) && !identical(rownames(object@junction_counts), cells)) {
      errors <- c(errors,
        "Row names of 'junction_counts' must match cell barcodes in the Seurat object.")
    }
  }

  # If the "psi" assay exists in Seurat, its columns must match cell barcodes
  if (!is.null(object@seurat) && inherits(object@seurat, "Seurat")) {
    psi_assay <- .get_assay_safe(object@seurat, "psi")
    if (!is.null(psi_assay)) {
      psi_cells <- colnames(psi_assay)
      if (!is.null(cells) && !identical(psi_cells, cells)) {
        errors <- c(errors,
          "Cell barcodes of the 'psi' assay must match those in the Seurat object.")
      }
      # Event IDs in "psi" assay must match event_data if both are present
      if (nrow(object@event_data) > 0) {
        psi_features <- rownames(psi_assay)
        if (!identical(sort(psi_features), sort(object@event_data$event_id))) {
          errors <- c(errors,
            "Features of the 'psi' assay must match 'event_id' in event_data.")
        }
      }
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

# Safe accessor for Seurat assays: returns NULL instead of erroring when absent.
# SeuratObject v5 [[]] throws an error if the assay name is not found.
.get_assay_safe <- function(seu, name) {
  if (is.null(seu) || !inherits(seu, "Seurat")) return(NULL)
  if (!name %in% SeuratObject::Assays(seu)) return(NULL)
  seu[[name]]
}

.get_cells <- function(object) {
  if (is.null(object@seurat)) return(NULL)
  colnames(object@seurat)
}

.n_cells <- function(object) {
  cells <- .get_cells(object)
  if (is.null(cells)) 0L else length(cells)
}

.n_events <- function(object) {
  if (is.null(object@seurat)) return(0L)
  psi_assay <- .get_assay_safe(object@seurat, "psi")
  if (is.null(psi_assay)) return(0L)
  nrow(psi_assay)
}
