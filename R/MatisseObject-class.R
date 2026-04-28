#' The MatisseObject S4 class
#'
#' The central data structure for Matisse. It wraps a \code{\link[Seurat]{Seurat}}
#' object and adds isoform-resolved splicing layers. All per-cell data -- junction
#' counts, PSI values, transcript counts, and QC metrics -- live inside the
#' embedded Seurat object as named assays (\code{Assay5}) or cell metadata
#' (\code{meta.data}). Nothing is duplicated outside the Seurat object.
#'
#' Two input modes are supported, set automatically at construction:
#' \describe{
#'   \item{\code{"junction"}}{Short-read mode. Raw junction counts (from
#'     STARsolo) are stored as \code{Assay5("isoform")} (junctions x cells).
#'     Call \code{\link{CalculatePSI}} to compute PSI values.}
#'   \item{\code{"transcript"}}{Long-read mode. Transcript isoform counts
#'     (from Bagpiper, FLAMES, or LIQA) are stored as
#'     \code{Assay5("isoform")} (transcripts x cells).
#'     Call \code{\link{CalculatePSI}} to compute PSI values.}
#' }
#'
#' After \code{\link{CalculatePSI}}, PSI values are stored as
#' \code{Assay5("psi")} (splice events x cells).
#'
#' @slot seurat A \code{Seurat} object. Contains all per-cell data: gene
#'   expression, splice assays (\code{"isoform"}, \code{"psi"}), cell metadata
#'   (QC metrics, cluster labels), and dimensionality reductions.
#' @slot input.mode Character. \code{"junction"} for short-read (STARsolo
#'   junction counts) objects; \code{"transcript"} for long-read (Bagpiper /
#'   FLAMES / LIQA transcript counts) objects.
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
    seurat       = "ANY",        # Seurat object
    input.mode   = "character",  # "junction" | "transcript"
    version      = "character",
    misc         = "list"
  ),
  prototype = list(
    seurat       = NULL,
    input.mode   = "junction",
    version      = as.character(utils::packageVersion("Matisse")),
    misc         = list()
  )
)

# ---------------------------------------------------------------------------
# Validity
# ---------------------------------------------------------------------------

setValidity("MatisseObject", function(object) {
  errors <- character()

  # input.mode must be one of the two supported values
  if (!object@input.mode %in% c("junction", "transcript")) {
    errors <- c(errors,
      "'input.mode' must be \"junction\" or \"transcript\".")
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
      # All event annotation event_ids must exist as features in the PSI assay
      ev <- .get_event_data_internal(object)
      if (!is.null(ev) && nrow(ev) > 0) {
        psi_features <- rownames(psi_assay)
        missing_ids  <- setdiff(ev$event_id, psi_features)
        if (length(missing_ids) > 0) {
          errors <- c(errors,
            paste0("Some event IDs in event annotation are missing from the ",
                   "'psi' assay: ",
                   paste(head(missing_ids, 3), collapse = ", ")))
        }
      }
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

# Internal helper: read event annotation from Misc() staging area.
# Returns a data.frame (or NULL if not present).
# Note: Misc<-.Seurat calls c() on list values, converting a data.frame
# to a named list. Reconstruct with as.data.frame().
.get_event_data_internal <- function(object) {
  if (is.null(object@seurat)) return(NULL)
  val <- SeuratObject::Misc(object@seurat, slot = "matisse_event_data")
  if (is.null(val)) return(NULL)
  as.data.frame(val, stringsAsFactors = FALSE)
}

.n_events <- function(object) {
  if (is.null(object@seurat)) return(0L)
  psi_assay <- .get_assay_safe(object@seurat, "psi")
  if (is.null(psi_assay)) return(0L)
  # event_data acts as a filter on the PSI assay (assay may hold a superset).
  # If event_data is present (even empty), honour it; otherwise count PSI rows.
  ev <- .get_event_data_internal(object)
  if (!is.null(ev)) {
    return(length(intersect(ev$event_id, rownames(psi_assay))))
  }
  nrow(psi_assay)
}

# Returns the number of raw input features stored in the "isoform" Assay5.
# Meaningful in both modes (junctions in junction mode, transcripts in
# transcript mode).
.n_junctions <- function(object) {
  if (object@input.mode != "junction") return(0L)
  iso_assay <- .get_assay_safe(object@seurat, "isoform")
  if (is.null(iso_assay)) return(0L)
  nrow(iso_assay)   # rows = features (Seurat's features x cells convention)
}
