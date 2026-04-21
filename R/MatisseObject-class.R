#' The MatisseObject S4 class
#'
#' The central data structure for Matisse. It wraps a \code{\link[Seurat]{Seurat}}
#' object and adds isoform-resolved splicing layers. All per-cell data — junction
#' counts, PSI values, transcript counts, and QC metrics — live inside the
#' embedded Seurat object as named assays (\code{Assay5}) or cell metadata
#' (\code{meta.data}). Nothing is duplicated outside the Seurat object.
#'
#' Two operating modes are supported, set automatically at construction:
#' \describe{
#'   \item{\code{"junction"}}{Short-read mode. Raw junction counts are stored
#'     as \code{Assay5("junction")} (junctions × cells). PSI is computed later
#'     by \code{\link{CalculatePSI}} and stored as \code{Assay5("psi")}.}
#'   \item{\code{"event"}}{Long-read mode. Transcript counts (e.g. from
#'     Bagpiper or FLAMES) are stored as \code{Assay5("transcript")}. PSI is
#'     computed at construction time from SUPPA2 \code{.ioe} event definitions
#'     and stored as \code{Assay5("psi")}.}
#' }
#'
#' @slot seurat A \code{Seurat} object. Contains all per-cell data: gene
#'   expression, splice assays (\code{"junction"}, \code{"transcript"},
#'   \code{"psi"}), cell metadata (QC metrics, cluster labels), and
#'   dimensionality reductions.
#' @slot event_data A \code{data.frame} with one row per splice event. Required
#'   columns: \code{event_id}, \code{gene_id}, \code{chr}, \code{strand},
#'   \code{event_type}, \code{inclusion_junctions}, \code{exclusion_junctions}.
#' @slot junction_data A \code{data.frame} with one row per junction. Required
#'   columns: \code{junction_id}, \code{chr}, \code{start}, \code{end},
#'   \code{strand}, \code{gene_id}.
#' @slot mode Character. \code{"junction"} for short-read objects;
#'   \code{"event"} for long-read objects.
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
    seurat        = "ANY",        # Seurat object
    event_data    = "data.frame",
    junction_data = "data.frame",
    mode          = "character",  # "junction" | "event"
    version       = "character",
    misc          = "list"
  ),
  prototype = list(
    seurat        = NULL,
    event_data    = data.frame(),
    junction_data = data.frame(),
    mode          = "junction",
    version       = as.character(utils::packageVersion("Matisse")),
    misc          = list()
  )
)

# ---------------------------------------------------------------------------
# Validity
# ---------------------------------------------------------------------------

setValidity("MatisseObject", function(object) {
  errors <- character()

  # mode must be one of the two supported values
  if (!object@mode %in% c("junction", "event")) {
    errors <- c(errors, "'mode' must be \"junction\" or \"event\".")
  }

  if (!is.null(object@seurat)) {
    if (!inherits(object@seurat, "Seurat")) {
      errors <- c(errors, "'seurat' slot must be a Seurat object or NULL.")
    }
  }

  cells <- .get_cells(object)

  # If the "psi" assay exists in Seurat, its columns must match cell barcodes
  if (!is.null(object@seurat) && inherits(object@seurat, "Seurat")) {
    psi_assay <- .get_assay_safe(object@seurat, "psi")
    if (!is.null(psi_assay)) {
      psi_cells <- colnames(psi_assay)
      if (!is.null(cells) && !identical(psi_cells, cells)) {
        errors <- c(errors,
          "Cell barcodes of the 'psi' assay must match those in the Seurat object.")
      }
      # All event_data event_ids must exist as features in the PSI assay
      # (the assay may hold a superset — e.g. after event-level subsetting)
      if (nrow(object@event_data) > 0) {
        psi_features <- rownames(psi_assay)
        missing_ids  <- setdiff(object@event_data$event_id, psi_features)
        if (length(missing_ids) > 0) {
          errors <- c(errors,
            paste0("Some event_data event_ids are missing from the 'psi' assay: ",
                   paste(head(missing_ids, 3), collapse = ", ")))
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
  # event_data is authoritative for active events (assay may hold a superset)
  nrow(object@event_data)
}

# Returns the number of junctions stored in the "junction" Assay5.
# Only meaningful in junction mode; returns 0L in event mode.
.n_junctions <- function(object) {
  if (object@mode != "junction") return(0L)
  jxn_assay <- .get_assay_safe(object@seurat, "junction")
  if (is.null(jxn_assay)) return(0L)
  nrow(jxn_assay)   # rows = junctions (Seurat's features × cells convention)
}
