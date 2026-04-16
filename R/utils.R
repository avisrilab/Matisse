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
#' via \code{merge()} (dispatches to Seurat's merge method).
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

  # Merge Seurat objects via base generic (dispatches to Seurat's merge.Seurat)
  merged_seurat <- merge(
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

# ---------------------------------------------------------------------------
# Shared PSI computation helpers
# ---------------------------------------------------------------------------

# Build a sparse (n_features x n_events) indicator matrix.
# id_lists   : list of length n_events; each element is a character vector of
#              feature IDs (transcript IDs or junction IDs) for that event.
# id_universe: character vector of all feature IDs (row universe = rownames of
#              the count matrix).
# Returns a dgCMatrix with 1s where feature i belongs to event j.
.build_indicator_matrix <- function(id_lists, id_universe) {
  n_ev  <- length(id_lists)
  n_ids <- length(id_universe)
  if (n_ev == 0L || n_ids == 0L) {
    return(Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = c(n_ids, n_ev), repr = "C"
    ))
  }
  id_vec  <- unlist(id_lists, use.names = FALSE)
  ev_idx  <- rep(seq_len(n_ev), lengths(id_lists))
  id_idx  <- match(id_vec, id_universe)
  keep    <- !is.na(id_idx)
  if (!any(keep)) {
    return(Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = c(n_ids, n_ev), repr = "C"
    ))
  }
  Matrix::sparseMatrix(
    i    = id_idx[keep],
    j    = ev_idx[keep],
    x    = 1.0,
    dims = c(n_ids, n_ev),
    repr = "C"
  )
}

# Compute a sparse PSI matrix from sparse inclusion and exclusion count matrices.
# inc_mat, exc_mat : dgCMatrix (cells x events)
# min_coverage     : integer
# Returns a dgCMatrix (cells x events) where only covered entries (total >=
# min_coverage) are stored.  Absent entries = not enough data (equivalent to NA).
# Use .psi_to_dense_na() when a dense NA-filled matrix is needed.
.psi_from_sparse_counts <- function(inc_mat, exc_mat, min_coverage) {
  total_mat <- inc_mat + exc_mat
  total_T   <- as(total_mat, "TsparseMatrix")   # triplet form
  nz_i      <- total_T@i                    # 0-based row
  nz_j      <- total_T@j                    # 0-based col
  nz_total  <- total_T@x
  n_cells   <- nrow(inc_mat)
  n_events  <- ncol(inc_mat)

  covered <- nz_total >= min_coverage

  # Look up inclusion counts at covered positions using column-major linear index.
  inc_T   <- as(inc_mat, "TsparseMatrix")
  inc_lin <- inc_T@j * as.double(n_cells) + inc_T@i   # 0-based col-major key

  psi_x <- rep(NA_real_, length(nz_total))             # default NA (low coverage)
  if (any(covered)) {
    cov_lin          <- nz_j[covered] * as.double(n_cells) + nz_i[covered]
    m                <- match(cov_lin, inc_lin)
    inc_cov          <- ifelse(is.na(m), 0.0, inc_T@x[m])
    psi_x[covered]   <- inc_cov / nz_total[covered]
  }
  # All non-zero-total entries are stored: covered → PSI value, not covered → NA.
  # Zero-total entries are absent (not stored); accessing them returns 0.
  Matrix::sparseMatrix(
    i    = nz_i + 1L,
    j    = nz_j + 1L,
    x    = psi_x,
    dims = c(n_cells, n_events),
    dimnames = dimnames(inc_mat),
    repr = "C"
  )
}

# Convert a sparse PSI matrix (only covered entries stored) to a dense matrix
# with NA for absent (uncovered) positions.  Used by SummarizePSI and any
# function that needs na.rm-aware statistics over cells.
.psi_to_dense_na <- function(psi_sparse) {
  psi_csc <- as(psi_sparse, "dgCMatrix")
  n_r     <- nrow(psi_sparse)
  n_c     <- ncol(psi_sparse)
  mat     <- matrix(NA_real_, nrow = n_r, ncol = n_c,
                    dimnames = dimnames(psi_sparse))
  n_stored <- length(psi_csc@i)
  if (n_stored > 0L) {
    col_idx <- rep(seq_len(n_c), diff(psi_csc@p))
    mat[cbind(psi_csc@i + 1L, col_idx)] <- psi_csc@x
  }
  mat
}

# Count covered cells per event: stored entries that are not NA.
# (NA entries = low coverage but non-zero total; absent entries = zero total.)
.n_covered_per_event <- function(psi_sparse) {
  psi_csc <- as(psi_sparse, "dgCMatrix")
  not_na  <- !is.na(psi_csc@x)
  col_idx <- rep(seq_len(ncol(psi_sparse)), diff(psi_csc@p))
  tabulate(col_idx[not_na], nbins = ncol(psi_sparse))
}

# Count covered events per cell: row-wise counterpart of .n_covered_per_event.
.n_covered_per_cell <- function(psi_sparse) {
  psi_csc <- as(psi_sparse, "dgCMatrix")
  not_na  <- !is.na(psi_csc@x)
  tabulate(psi_csc@i[not_na] + 1L, nbins = nrow(psi_sparse))
}

# Sparse row means of PSI, ignoring NAs, without any dense conversion.
.psi_rowmeans_sparse <- function(psi_sparse) {
  psi_csc      <- as(psi_sparse, "dgCMatrix")
  psi_x        <- psi_csc@x
  psi_x[is.na(psi_x)] <- 0
  psi_csc@x    <- psi_x
  n_cov        <- .n_covered_per_cell(psi_sparse)
  row_sums     <- as.numeric(Matrix::rowSums(psi_csc))
  ifelse(n_cov > 0L, row_sums / n_cov, NA_real_)
}
