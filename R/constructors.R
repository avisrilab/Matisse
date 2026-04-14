#' @include MatisseObject-class.R
#' @include generics.R
NULL

#' Create a MatisseObject
#'
#' The primary constructor for \code{\linkS4class{MatisseObject}}. Combines a
#' \code{Seurat} object with optional isoform-resolved data layers.
#'
#' @param seurat A \code{Seurat} object. Required; provides cell barcodes,
#'   gene expression, and cell-level metadata.
#' @param junction_counts A sparse matrix (dgCMatrix, cells x junctions) of
#'   raw per-junction read counts. Row names must match
#'   \code{colnames(seurat)}. Default: \code{NULL}.
#' @param event_data A \code{data.frame} defining splice events. Required
#'   columns: \code{event_id}, \code{gene_id}, \code{chr}, \code{strand},
#'   \code{event_type}, \code{inclusion_junctions} (semicolon-separated
#'   junction IDs), \code{exclusion_junctions} (semicolon-separated junction
#'   IDs). Default: \code{NULL} (empty table).
#' @param junction_data A \code{data.frame} of junction annotations. Required
#'   columns: \code{junction_id}, \code{chr}, \code{start}, \code{end},
#'   \code{strand}, \code{gene_id}. Default: \code{NULL} (empty table).
#' @param verbose Logical. Print construction progress. Default: \code{TRUE}.
#'
#' @return A \code{\linkS4class{MatisseObject}}.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' counts <- matrix(rpois(200, 5), nrow = 20,
#'                  dimnames = list(paste0("Gene", 1:20),
#'                                  paste0("Cell", 1:10)))
#' seu <- CreateSeuratObject(counts)
#' obj <- CreateMatisseObject(seurat = seu)
#' obj
#' }
#'
#' @export
CreateMatisseObject <- function(
    seurat,
    junction_counts = NULL,
    event_data      = NULL,
    junction_data   = NULL,
    verbose         = TRUE
) {
  # --- validate seurat -------------------------------------------------------
  if (!inherits(seurat, "Seurat")) {
    rlang::abort("`seurat` must be a Seurat object.")
  }
  cells <- colnames(seurat)

  if (verbose) {
    cli::cli_alert_info(
      "Creating MatisseObject: {length(cells)} cells, \\
       {nrow(seurat)} gene features.")
  }

  # --- validate junction_counts ----------------------------------------------
  if (!is.null(junction_counts)) {
    junction_counts <- .validate_cell_matrix(
      junction_counts, cells, "junction_counts")
  }

  # --- validate event_data ---------------------------------------------------
  if (is.null(event_data)) {
    event_data <- data.frame()
  } else {
    required <- c("event_id", "gene_id", "chr", "strand",
                  "event_type", "inclusion_junctions", "exclusion_junctions")
    .check_required_columns(event_data, required, "event_data")
    event_data <- as.data.frame(event_data)
  }

  # --- validate junction_data ------------------------------------------------
  if (is.null(junction_data)) {
    junction_data <- data.frame()
  } else {
    required <- c("junction_id", "chr", "start", "end", "strand", "gene_id")
    .check_required_columns(junction_data, required, "junction_data")
    junction_data <- as.data.frame(junction_data)
  }

  # --- build isoform_metadata skeleton ---------------------------------------
  isoform_metadata <- data.frame(row.names = cells)

  # --- assemble ---------------------------------------------------------------
  version_str <- tryCatch(
    as.character(utils::packageVersion("Matisse")),
    error = function(e) "development"
  )

  obj <- methods::new(
    "MatisseObject",
    seurat           = seurat,
    psi              = NULL,
    inclusion_counts = NULL,
    exclusion_counts = NULL,
    junction_counts  = junction_counts,
    event_data       = event_data,
    junction_data    = junction_data,
    isoform_metadata = isoform_metadata,
    version          = version_str,
    misc             = list()
  )

  if (verbose) {
    cli::cli_alert_success("MatisseObject created successfully.")
  }

  methods::validObject(obj)
  obj
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Coerce to dgCMatrix and check row names align with cell barcodes
.validate_cell_matrix <- function(mat, cells, name) {
  if (!inherits(mat, "Matrix")) {
    mat <- Matrix::Matrix(as.matrix(mat), sparse = TRUE)
  }
  if (is.null(rownames(mat))) {
    rlang::abort(paste0("`", name, "` must have row names matching cell barcodes."))
  }
  if (!identical(sort(rownames(mat)), sort(cells))) {
    # Allow partial overlap with a warning; reorder to match Seurat order
    common <- intersect(rownames(mat), cells)
    if (length(common) == 0) {
      rlang::abort(paste0(
        "No row names of `", name, "` match cell barcodes in `seurat`."))
    }
    if (length(common) < length(cells)) {
      rlang::warn(paste0(
        "`", name, "` covers ", length(common), "/", length(cells),
        " cells. Missing cells will have all-zero rows."))
      # Pad with zeros for missing cells
      missing_cells <- setdiff(cells, rownames(mat))
      pad <- Matrix::sparseMatrix(
        i = integer(0), j = integer(0),
        dims = c(length(missing_cells), ncol(mat)),
        dimnames = list(missing_cells, colnames(mat))
      )
      mat <- rbind(mat[common, ], pad)
    }
    mat <- mat[cells, , drop = FALSE]
  }
  methods::as(mat, "CsparseMatrix")
}

# Abort if a data.frame is missing required column names
.check_required_columns <- function(df, required, name) {
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    rlang::abort(paste0(
      "`", name, "` is missing required columns: ",
      paste(missing, collapse = ", ")))
  }
  invisible(TRUE)
}
