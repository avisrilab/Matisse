#' @include MatisseObject-class.R
NULL

# ---------------------------------------------------------------------------
# Exported utilities
# ---------------------------------------------------------------------------

#' Merge two MatisseObjects by cells
#'
#' Concatenates two \code{MatisseObject}s that share the same set of splice
#' events. The embedded Seurat objects are merged via \code{merge()} (which
#' dispatches to Seurat's merge method, combining all assays including
#' \code{"isoform"} and \code{"psi"}). Both objects must have the same
#' \code{input.mode}.
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
  if (x@input.mode != y@input.mode) {
    rlang::abort(paste0(
      "Cannot merge MatisseObjects with different input.mode values: ",
      "x is '", x@input.mode, "', y is '", y@input.mode, "'."))
  }
  stopifnot(length(add_cell_ids) == 2L,
            is.character(add_cell_ids),
            all(nzchar(add_cell_ids)))

  # Check event compatibility via the "psi" assay features
  psi_x <- GetPSI(x)
  psi_y <- GetPSI(y)
  if (!is.null(psi_x) && !is.null(psi_y)) {
    if (!identical(colnames(psi_x), colnames(psi_y))) {
      rlang::warn(
        "PSI event names differ between objects. Only shared events will be kept.")
      shared <- intersect(colnames(psi_x), colnames(psi_y))
      x <- x[, shared]
      y <- y[, shared]
    }
  }

  # Merge Seurat objects -- handles all assays (isoform, psi, RNA) and cell
  # metadata automatically. JoinLayers consolidates per-sample layers back into
  # unified layers (SeuratObject v5 splits them on merge). The PSI assay's
  # meta.features (event annotation) is also carried through.
  merged_seurat <- merge(
    x@seurat, y@seurat,
    add.cell.ids = add_cell_ids
  )
  for (assay_name in SeuratObject::Assays(merged_seurat)) {
    merged_seurat <- SeuratObject::JoinLayers(merged_seurat, assay = assay_name)
  }

  # @misc holds provenance/staging fields. Use modifyList so x wins per-key
  # overall (matching the docstring), not just for the two keys we'd
  # otherwise touch by hand.
  merged_misc <- utils::modifyList(y@misc, x@misc)

  obj <- methods::new(
    "MatisseObject",
    seurat     = merged_seurat,
    input.mode = x@input.mode,
    version    = x@version,
    misc       = merged_misc
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

# ---------------------------------------------------------------------------
# Junction-name coordinate parser
# ---------------------------------------------------------------------------

# Try to parse genomic coordinates out of junction IDs. Tools encode them
# differently — STARsolo's "chr1-12345-67890-+", SUPPA-style
# "chr1:12345-67890:+", underscore-separated "chr1_12345_67890_+", etc.
# Returns a data.frame keyed by `junction_id` with columns chr/start/end/strand.
# Names that don't match any pattern get NA coordinates and produce a warning;
# everything downstream still works for raw counts, only sashimi degrades.
.parse_junction_names <- function(junction_ids) {
  n <- length(junction_ids)
  out <- data.frame(
    junction_id = junction_ids,
    chr         = NA_character_,
    start       = NA_integer_,
    end         = NA_integer_,
    strand      = NA_character_,
    stringsAsFactors = FALSE
  )
  if (n == 0L) return(out)

  patterns <- c(
    # chr1-12345-67890-+   (STARsolo SJ default)
    "^([^-]+)-(\\d+)-(\\d+)-([-+*.])$",
    # chr1:12345-67890:+   (SUPPA / common bedlike)
    "^([^:]+):(\\d+)-(\\d+):([-+*.])$",
    # chr1_12345_67890_+   (underscore-separated)
    "^([^_]+)_(\\d+)_(\\d+)_([-+*.])$"
  )

  unmatched <- rep(TRUE, n)
  for (pat in patterns) {
    if (!any(unmatched)) break
    candidates <- junction_ids[unmatched]
    matches    <- regmatches(candidates, regexec(pat, candidates))
    hit_local  <- vapply(matches, length, integer(1)) == 5L
    if (!any(hit_local)) next
    abs_idx    <- which(unmatched)[hit_local]
    parts      <- matches[hit_local]
    out$chr[abs_idx]    <- vapply(parts, function(x) x[2L], character(1))
    out$start[abs_idx]  <- as.integer(vapply(parts, function(x) x[3L], character(1)))
    out$end[abs_idx]    <- as.integer(vapply(parts, function(x) x[4L], character(1)))
    out$strand[abs_idx] <- vapply(parts, function(x) x[5L], character(1))
    unmatched[abs_idx]  <- FALSE
  }

  if (any(unmatched)) {
    n_bad <- sum(unmatched)
    examples <- paste(utils::head(junction_ids[unmatched], 3L), collapse = ", ")
    rlang::warn(paste0(
      "Could not parse coordinates from ", n_bad, " junction ID(s); ",
      "they will have NA coords (e.g. ", examples, "). ",
      "Junction-mode PlotSashimi may degrade for these events. ",
      "Supported formats: chr-start-end-strand, chr:start-end:strand, ",
      "chr_start_end_strand."))
  }
  out
}

# ---------------------------------------------------------------------------
# Shared PSI computation helpers
# ---------------------------------------------------------------------------

.build_indicator_matrix <- function(id_lists, id_universe) {
  n_ev  <- length(id_lists)
  n_ids <- length(id_universe)
  if (n_ev == 0L || n_ids == 0L) {
    return(Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = c(n_ids, n_ev), repr = "C"
    ))
  }
  id_vec <- unlist(id_lists, use.names = FALSE)
  ev_idx <- rep(seq_len(n_ev), lengths(id_lists))
  id_idx <- match(id_vec, id_universe)
  keep   <- !is.na(id_idx)
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

.psi_from_sparse_counts <- function(inc_mat, exc_mat, min_coverage,
                                     n_inc = NULL, n_exc = NULL) {
  total_mat <- inc_mat + exc_mat
  total_T   <- as(total_mat, "TsparseMatrix")
  nz_i      <- total_T@i
  nz_j      <- total_T@j
  nz_total  <- total_T@x
  n_cells   <- nrow(inc_mat)
  n_events  <- ncol(inc_mat)

  # Coverage check uses RAW total event-supporting reads regardless of
  # whether the PSI ratio is normalized below — keeps min_coverage's
  # user-facing meaning ("at least N reads at this event in this cell")
  # stable across the two paths.
  covered <- nz_total >= min_coverage

  inc_T   <- as(inc_mat, "TsparseMatrix")
  inc_lin <- inc_T@j * as.double(n_cells) + inc_T@i

  psi_x <- rep(NA_real_, length(nz_total))
  if (any(covered)) {
    cov_lin        <- nz_j[covered] * as.double(n_cells) + nz_i[covered]
    m              <- match(cov_lin, inc_lin)
    inc_cov        <- ifelse(is.na(m), 0.0, inc_T@x[m])

    if (!is.null(n_inc) && !is.null(n_exc)) {
      # SUPPA2 / rMATS-style per-junction normalization. Without this,
      # events with asymmetric junction counts (SE has 2 inclusion vs.
      # 1 exclusion) over-weight the side with more junctions because
      # the same molecule produces reads at multiple junctions on that
      # side. Dividing each role's read sum by its junction count
      # recovers the molecular fraction.
      exc_cov        <- nz_total[covered] - inc_cov
      event_idx      <- nz_j[covered] + 1L
      inc_norm       <- inc_cov / n_inc[event_idx]
      exc_norm       <- exc_cov / n_exc[event_idx]
      psi_x[covered] <- inc_norm / (inc_norm + exc_norm)
    } else {
      psi_x[covered] <- inc_cov / nz_total[covered]
    }
  }
  Matrix::sparseMatrix(
    i    = nz_i + 1L,
    j    = nz_j + 1L,
    x    = psi_x,
    dims = c(n_cells, n_events),
    dimnames = dimnames(inc_mat),
    repr = "C"
  )
}

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

.n_covered_per_event <- function(psi_sparse) {
  psi_csc <- as(psi_sparse, "dgCMatrix")
  not_na  <- !is.na(psi_csc@x)
  col_idx <- rep(seq_len(ncol(psi_sparse)), diff(psi_csc@p))
  tabulate(col_idx[not_na], nbins = ncol(psi_sparse))
}

.n_covered_per_cell <- function(psi_sparse) {
  psi_csc <- as(psi_sparse, "dgCMatrix")
  not_na  <- !is.na(psi_csc@x)
  tabulate(psi_csc@i[not_na] + 1L, nbins = nrow(psi_sparse))
}
