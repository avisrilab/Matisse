#' @include MatisseObject-class.R
#' @include generics.R
NULL

# ---------------------------------------------------------------------------
# CalculatePSI -- MatisseObject method
# ---------------------------------------------------------------------------

#' Calculate PSI matrix from junction or transcript counts
#'
#' Computes a Percent Spliced In (PSI) matrix for all splice events and stores
#' it in the \code{"psi"} assay. Works in both junction mode and transcript
#' mode. \code{\link{CreateMatisseObject}} calls this automatically by
#' default; you only call it directly to recompute with different parameters
#' (e.g. a different \code{min_coverage}) or after constructing with
#' \code{defer_psi = TRUE}.
#'
#' For each cell \eqn{c} and event \eqn{e}:
#'
#' \deqn{PSI_{c,e} = \frac{\sum \text{inclusion reads}}
#'                        {\sum \text{inclusion reads} +
#'                         \sum \text{exclusion reads}}}
#'
#' Results are stored inside the embedded Seurat object as
#' \code{Assay5("psi")}, with:
#' \itemize{
#'   \item \code{"data"} layer: PSI values in \eqn{[0,1]} (events x cells).
#'   \item \code{"counts"} layer: inclusion read counts (events x cells).
#'   \item \code{"exclusion"} layer: exclusion read counts (events x cells).
#' }
#' Entries where total coverage falls below \code{min_coverage} are set to
#' \code{NA} in the \code{"data"} layer.
#'
#' \code{nPercent_isoform} is also written to cell metadata: the percentage of
#' splice events with a non-\code{NA} PSI value in each cell.
#'
#' @param object A \code{\linkS4class{MatisseObject}}, or a sparse matrix
#'   (cells x junctions) when computing PSI outside the object.
#' @param events A \code{data.frame} with columns \code{event_id},
#'   \code{inclusion_features}, and \code{exclusion_features}. Defaults to
#'   the event annotation staged at construction time
#'   (\code{object@misc[["event_data"]]}, populated by
#'   \code{\link{CreateMatisseObject}}).
#' @param min_coverage Integer. Minimum total reads per cell per event to
#'   report a PSI value. Default: \code{5}.
#' @param verbose Logical. Print progress. Default: \code{TRUE}.
#'
#' @return
#' * \code{MatisseObject}: the input object with the \code{"psi"} assay
#'   populated inside the embedded Seurat object.
#' * matrix: a dense matrix (cells x events) of PSI values.
#'
#' @seealso \code{\link{FilterCells}}, \code{\link{FilterEvents}},
#'   \code{\link{PlotHeatmap}}
#'
#' @rdname CalculatePSI
#' @export
setMethod("CalculatePSI", "MatisseObject",
          function(object, events = NULL, min_coverage = 5L,
                   verbose = TRUE, ...) {

  iso_assay <- .get_assay_safe(object@seurat, "isoform")

  if (object@input.mode == "junction") {
    # --- Junction mode: aggregate junction reads to splice events ------------
    if (is.null(iso_assay)) {
      rlang::abort(paste0(
        "No junction assay found. ",
        "Provide junction_counts via CreateMatisseObject()."))
    }
    # Assay5 stores junctions x cells; Matisse internal convention is
    # cells x junctions, so transpose.
    jxn_counts <- Matrix::t(.get_assay_layer(iso_assay, "counts"))
  }

  if (is.null(events)) {
    # First read from the @misc staging area (populated at construction).
    # If empty, fall back to the PSI assay's meta.features (CalculatePSI was
    # called previously and migrated event_data into the assay).
    events <- object@misc[["event_data"]]
    if (is.null(events) || nrow(events) == 0) {
      psi_assay <- .get_assay_safe(object@seurat, "psi")
      if (!is.null(psi_assay)) {
        mf <- psi_assay[[]]
        if (!is.null(mf) && nrow(mf) > 0L) {
          mf$event_id <- rownames(mf)
          events <- mf
        }
      }
    }
  }
  if (is.null(events) || nrow(events) == 0) {
    rlang::abort(paste0(
      "No splice events defined. Provide events via CreateMatisseObject() ",
      "or the `events` argument."))
  }

  if (object@input.mode == "junction") {
    result <- .calculate_psi_matrix(
      jxn_counts   = jxn_counts,
      events       = events,
      min_coverage = min_coverage,
      verbose      = verbose
    )
  } else {
    # --- Transcript mode: aggregate transcript counts to splice events -------
    if (is.null(iso_assay)) {
      rlang::abort(paste0(
        "No 'isoform' assay found. ",
        "Provide transcript_counts via CreateMatisseObject()."))
    }
    tx_counts <- .get_assay_layer(iso_assay, "counts")
    if (verbose) {
      cli::cli_alert_info(paste0(
        "Calculating PSI for {nrow(events)} events across ",
        "{ncol(tx_counts)} cells..."))
    }
    result <- .aggregate_transcript_counts(
      tx_counts    = tx_counts,
      events       = events,
      min_coverage = min_coverage,
      cells        = colnames(tx_counts)
    )
  }

  # Store PSI, inclusion, and exclusion in an Assay5 named "psi"
  psi_result <- .create_psi_assay(
    psi_mat = result$psi,
    inc_mat = result$inclusion,
    exc_mat = result$exclusion
  )
  object@seurat[["psi"]] <- psi_result$assay

  # Migrate event annotation into the PSI assay's feature metadata. After
  # CalculatePSI returns, event_data lives in seurat[["psi"]][[]] and is
  # the single source of truth — kept in sync with rownames automatically by
  # Seurat under any subset/merge. The transient @misc[["event_data"]]
  # staging area populated by CreateMatisseObject is cleared.
  stored_names <- psi_result$feature_names
  ev <- as.data.frame(events, stringsAsFactors = FALSE)
  ev$event_id  <- stored_names
  rownames(ev) <- stored_names
  # Seurat's [[<-.Assay5 requires rownames matching the assay's feature names;
  # sanitization above guarantees this.
  object@seurat[["psi"]][[]] <- ev
  object@misc[["event_data"]] <- NULL

  # Write nPercent_isoform to meta.data using the sparse-aware coverage
  # counter — avoids the dense (cells x events) materialisation that
  # .psi_to_dense_na would otherwise force.
  psi_cx  <- GetPSI(object)         # cells x events (sparse)
  n_total <- ncol(psi_cx)
  if (n_total > 0L) {
    n_cov_per_cell <- .n_covered_per_cell(psi_cx)
    meta_df <- data.frame(
      nPercent_isoform = round(100 * n_cov_per_cell / n_total, 2),
      row.names        = rownames(psi_cx)
    )
    object@seurat <- SeuratObject::AddMetaData(object@seurat, meta_df)
  }

  if (verbose) {
    psi_sp <- result$psi
    pct <- round(100 * sum(.n_covered_per_event(psi_sp)) /
                   (as.double(nrow(psi_sp)) * ncol(psi_sp)), 1)
    cli::cli_alert_success(
      "PSI calculated for {ncol(psi_sp)} events; {pct}% entries covered.")
  }

  methods::validObject(object)
  object
})

# CalculatePSI for a raw matrix input
#' @rdname CalculatePSI
setMethod("CalculatePSI", "ANY",
          function(object, events, min_coverage = 5L,
                   verbose = TRUE, ...) {
  if (!inherits(object, "Matrix") && !is.matrix(object)) {
    rlang::abort("`object` must be a MatisseObject or a matrix.")
  }
  if (is.null(events) || nrow(events) == 0) {
    rlang::abort("`events` is required when `object` is a matrix.")
  }

  result <- .calculate_psi_matrix(
    jxn_counts   = object,
    events       = events,
    min_coverage = min_coverage,
    verbose      = verbose
  )
  .psi_to_dense_na(result$psi)
})

# ---------------------------------------------------------------------------
# Core PSI computation (internal)
# ---------------------------------------------------------------------------

.calculate_psi_matrix <- function(jxn_counts, events,
                                   min_coverage, verbose) {
  .check_required_columns(
    events,
    c("event_id", "inclusion_features", "exclusion_features"),
    "events"
  )

  jxn_names <- colnames(jxn_counts)
  n_cells   <- nrow(jxn_counts)
  n_events  <- nrow(events)
  cells     <- rownames(jxn_counts)

  if (verbose) {
    cli::cli_alert_info(
      "Calculating PSI for {n_events} events across {n_cells} cells...")
  }

  inc_lists <- lapply(strsplit(events$inclusion_features, ";", fixed = TRUE),
                      trimws)
  exc_lists <- lapply(strsplit(events$exclusion_features, ";", fixed = TRUE),
                      trimws)

  A_inc <- .build_indicator_matrix(inc_lists, jxn_names)
  A_exc <- .build_indicator_matrix(exc_lists, jxn_names)
  colnames(A_inc) <- colnames(A_exc) <- events$event_id

  inc_mat <- jxn_counts %*% A_inc
  exc_mat <- jxn_counts %*% A_exc
  dimnames(inc_mat) <- dimnames(exc_mat) <- list(cells, events$event_id)

  psi_mat <- .psi_from_sparse_counts(inc_mat, exc_mat, min_coverage)

  list(
    psi       = psi_mat,
    inclusion = Matrix::Matrix(round(inc_mat), sparse = TRUE),
    exclusion = Matrix::Matrix(round(exc_mat), sparse = TRUE)
  )
}

# ---------------------------------------------------------------------------
# SummarizePSI
# ---------------------------------------------------------------------------

#' Summarize PSI distribution across cells for each event
#'
#' Returns a summary table with per-event PSI statistics across all (or a
#' subset of) cells. Call this after \code{\link{CalculatePSI}}.
#'
#' @param object A \code{MatisseObject} with a \code{"psi"} assay.
#' @param cells Optional character vector of cell barcodes to subset.
#'
#' @return A \code{data.frame} with one row per event and columns:
#'   \code{event_id}, \code{mean_psi}, \code{median_psi}, \code{sd_psi},
#'   \code{n_cells_covered}.
#'
#' @export
SummarizePSI <- function(object, cells = NULL) {
  if (!inherits(object, "MatisseObject")) {
    rlang::abort("`object` must be a MatisseObject.")
  }
  psi_cx <- GetPSI(object)  # cells x events
  if (is.null(psi_cx)) {
    rlang::abort("PSI matrix is empty. Run CalculatePSI() first.")
  }

  if (!is.null(cells)) psi_cx <- psi_cx[cells, , drop = FALSE]

  psi_csc <- as(psi_cx, "dgCMatrix")       # cells x events, sparse
  x       <- psi_csc@x
  not_na  <- !is.na(x)
  col_idx <- rep(seq_len(ncol(psi_csc)), diff(psi_csc@p))

  # Coverage (number of non-NA entries per event column)
  n_cov <- tabulate(col_idx[not_na], nbins = ncol(psi_csc))

  # Mean: colSums(non-NA values) / n_cov -- no dense matrix
  x_zeroed        <- x; x_zeroed[!not_na] <- 0
  psi_for_sum     <- psi_csc; psi_for_sum@x <- x_zeroed
  col_sums        <- as.numeric(Matrix::colSums(psi_for_sum))
  mean_psi        <- ifelse(n_cov > 0L, col_sums / n_cov, NA_real_)

  # SD: Var = (sum(x^2) - sum(x)^2 / n) / (n-1), SD = sqrt(Var)
  psi_sq          <- psi_for_sum; psi_sq@x <- x_zeroed^2
  col_sum_sq      <- as.numeric(Matrix::colSums(psi_sq))
  sd_psi          <- ifelse(
    n_cov > 1L,
    sqrt(pmax(0, (col_sum_sq - col_sums^2 / n_cov) / (n_cov - 1L))),
    NA_real_
  )

  # Median: tapply over the non-NA stored values -- avoids full dense alloc
  vals_non_na <- x[not_na]
  cols_non_na <- col_idx[not_na]
  median_psi  <- rep(NA_real_, ncol(psi_csc))
  if (length(vals_non_na) > 0L) {
    med_vals   <- tapply(vals_non_na, cols_non_na, stats::median)
    median_psi[as.integer(names(med_vals))] <- as.numeric(med_vals)
  }

  data.frame(
    event_id        = colnames(psi_csc),
    mean_psi        = mean_psi,
    median_psi      = median_psi,
    sd_psi          = sd_psi,
    n_cells_covered = n_cov,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
}
