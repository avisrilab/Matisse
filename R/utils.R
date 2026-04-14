#' @include MatisseObject-class.R
NULL

# ---------------------------------------------------------------------------
# Exported utilities
# ---------------------------------------------------------------------------

#' Build a minimal junction event annotation table
#'
#' Helper for quickly creating the \code{event_data} data.frame required by
#' \code{\link{CreateMatisseObject}} and \code{\link{CalculatePSI}} from a
#' table of per-junction metadata.
#'
#' Each row in \code{junctions} defines one junction. Events are then built
#' by pairing junctions from the same gene according to \code{event_type}:
#' \describe{
#'   \item{simple}{Each junction is treated as its own event. Inclusion =
#'     the junction itself; exclusion = all other junctions for the same gene.}
#' }
#'
#' This is a convenience function for exploratory use. For production
#' analyses, supply a fully annotated \code{event_data} table (e.g. from
#' rMATS or SUPPA2).
#'
#' @param junctions A \code{data.frame} with at minimum columns
#'   \code{junction_id} and \code{gene_id}.
#' @param event_type Character. Currently only \code{"simple"} is supported.
#'
#' @return A \code{data.frame} suitable for the \code{event_data} argument
#'   of \code{\link{CreateMatisseObject}}.
#'
#' @export
BuildSimpleEvents <- function(junctions, event_type = "simple") {
  .check_required_columns(junctions, c("junction_id", "gene_id"), "junctions")

  if (event_type != "simple") {
    rlang::abort("Only event_type = 'simple' is currently supported.")
  }

  genes <- unique(junctions$gene_id)
  rows  <- lapply(genes, function(g) {
    jxns <- junctions$junction_id[junctions$gene_id == g]
    lapply(jxns, function(focal) {
      other <- setdiff(jxns, focal)
      data.frame(
        event_id             = paste0(g, ":", focal),
        gene_id              = g,
        chr                  = .safe_col(junctions, focal, "chr"),
        strand               = .safe_col(junctions, focal, "strand"),
        event_type           = "simple",
        inclusion_junctions  = focal,
        exclusion_junctions  = paste(other, collapse = ";"),
        stringsAsFactors     = FALSE
      )
    })
  })

  do.call(rbind, unlist(rows, recursive = FALSE))
}

#' Merge two MatisseObjects by cells
#'
#' Concatenates two \code{MatisseObject}s that share the same set of splice
#' events. All matrices are row-bound; the embedded Seurat objects are merged
#' with \code{\link[SeuratObject]{merge}}.
#'
#' @param x A \code{MatisseObject}.
#' @param y A \code{MatisseObject}.
#' @param add_cell_ids Character vector of length 2. Prefixes appended to
#'   cell barcodes to avoid collisions. Default: \code{c("x", "y")}.
#' @param verbose Logical. Default: \code{TRUE}.
#'
#' @return A merged \code{MatisseObject}.
#'
#' @export
MergeMatisse <- function(x, y, add_cell_ids = c("x", "y"), verbose = TRUE) {
  if (!inherits(x, "MatisseObject") || !inherits(y, "MatisseObject")) {
    rlang::abort("Both `x` and `y` must be MatisseObjects.")
  }

  # Check event compatibility
  if (!is.null(x@psi) && !is.null(y@psi)) {
    if (!identical(colnames(x@psi), colnames(y@psi))) {
      rlang::warn(
        "PSI column names differ between objects. Only shared events will be kept.")
      shared <- intersect(colnames(x@psi), colnames(y@psi))
      x <- x[, shared]
      y <- y[, shared]
    }
  }

  # Merge Seurat objects
  merged_seurat <- SeuratObject::merge(
    x@seurat, y@seurat,
    add.cell.ids = add_cell_ids
  )

  rbind_sparse <- function(m1, m2, prefix1, prefix2) {
    if (is.null(m1) || is.null(m2)) return(NULL)
    rownames(m1) <- paste0(prefix1, "_", rownames(m1))
    rownames(m2) <- paste0(prefix2, "_", rownames(m2))
    rbind(m1, m2)
  }

  p1 <- add_cell_ids[1]; p2 <- add_cell_ids[2]

  obj <- methods::new(
    "MatisseObject",
    seurat           = merged_seurat,
    psi              = rbind_sparse(x@psi,              y@psi,              p1, p2),
    inclusion_counts = rbind_sparse(x@inclusion_counts, y@inclusion_counts, p1, p2),
    exclusion_counts = rbind_sparse(x@exclusion_counts, y@exclusion_counts, p1, p2),
    junction_counts  = rbind_sparse(x@junction_counts,  y@junction_counts,  p1, p2),
    event_data       = x@event_data,
    junction_data    = x@junction_data,
    isoform_metadata = data.frame(
      row.names = colnames(merged_seurat)),
    version          = x@version,
    misc             = c(x@misc, y@misc)
  )

  if (verbose) {
    cli::cli_alert_success(
      "Merged {.n_cells(x)} + {.n_cells(y)} = {.n_cells(obj)} cells.")
  }

  obj
}

# ---------------------------------------------------------------------------
# Internal utilities (not exported)
# ---------------------------------------------------------------------------

# Safely look up a column from junction_data for a specific junction_id
.safe_col <- function(df, id, col) {
  if (!col %in% colnames(df)) return(NA_character_)
  val <- df[[col]][df$junction_id == id]
  if (length(val) == 0) NA_character_ else as.character(val[1])
}

# Check required column names (re-exported from constructors.R for reuse)
# (defined in constructors.R; referenced here via the package namespace)
