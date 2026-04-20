#' @include MatisseObject-class.R
#' @include generics.R
NULL

# ---------------------------------------------------------------------------
# ComputeIsoformQC
# ---------------------------------------------------------------------------

#' Compute per-cell isoform quality-control metrics
#'
#' Calculates a panel of QC metrics from the junction count and PSI layers and
#' stores them in the \code{isoform_metadata} slot. Existing columns with the
#' same names are overwritten.
#'
#' Computed metrics:
#' \describe{
#'   \item{n_junctions_detected}{Number of junctions with at least one read.}
#'   \item{total_junction_reads}{Total junction read count across all junctions.}
#'   \item{n_events_covered}{Number of splice events with PSI \eqn{\geq}
#'     \code{min_coverage} (requires PSI to be calculated first).}
#'   \item{pct_events_covered}{Percentage of events with sufficient coverage.}
#'   \item{mean_psi}{Mean PSI across covered events.}
#' }
#'
#' @param object A \code{MatisseObject}.
#' @param min_coverage Integer. Minimum total reads to consider an event
#'   "covered" when computing \code{n_events_covered}. Default: \code{5}.
#' @param verbose Logical. Default: \code{TRUE}.
#'
#' @return The \code{MatisseObject} with QC columns added to
#'   \code{isoform_metadata}.
#'
#' @seealso \code{\link{FilterCells}}, \code{\link{PlotQCMetrics}}
#'
#' @rdname ComputeIsoformQC
#' @export
setMethod("ComputeIsoformQC", "MatisseObject",
          function(object, min_coverage = 5L, verbose = TRUE) {

  qc <- data.frame(row.names = .get_cells(object))

  # --- junction-level metrics ------------------------------------------------
  if (!is.null(object@junction_counts)) {
    jxn <- object@junction_counts
    qc$n_junctions_detected <- as.integer(Matrix::rowSums(jxn > 0))
    qc$total_junction_reads <- as.integer(Matrix::rowSums(jxn))
  } else {
    rlang::warn(
      "junction_counts is NULL; junction QC metrics will not be computed.")
  }

  # --- PSI-level metrics — retrieve from the "psi" assay --------------------
  psi_cx <- GetPSI(object)  # cells x events
  if (!is.null(psi_cx)) {
    n_events <- ncol(psi_cx)
    inc_cx   <- GetInclusionCounts(object)  # cells x events
    exc_cx   <- GetExclusionCounts(object)  # cells x events

    if (!is.null(inc_cx) && !is.null(exc_cx)) {
      total_csc   <- as(inc_cx + exc_cx, "dgCMatrix")
      total_csc@x <- as.numeric(total_csc@x >= min_coverage)
      covered_sp  <- Matrix::drop0(total_csc)
      n_cov       <- as.integer(Matrix::rowSums(covered_sp))
      qc$n_events_covered   <- n_cov
      qc$pct_events_covered <- round(100 * n_cov / n_events, 2)
    }

    qc$mean_psi <- .psi_rowmeans_sparse(psi_cx)
  } else {
    if (verbose) {
      rlang::warn(
        "PSI matrix is NULL; PSI-level QC metrics will not be computed. ",
        "Run CalculatePSI() first.")
    }
  }

  object <- AddIsoformMetadata(object, qc)

  if (verbose) {
    cli::cli_alert_success(
      "Computed {ncol(qc)} QC metrics for {nrow(qc)} cells.")
  }

  object
})

# ---------------------------------------------------------------------------
# FilterCells
# ---------------------------------------------------------------------------

#' Filter cells by isoform QC thresholds
#'
#' Removes cells that do not pass the specified thresholds on columns in the
#' \code{isoform_metadata} slot.
#'
#' @param object A \code{MatisseObject}.
#' @param min_junctions Integer. Minimum \code{n_junctions_detected}.
#'   Default: \code{NULL} (no filter).
#' @param max_junctions Integer. Maximum \code{n_junctions_detected}.
#'   Default: \code{NULL}.
#' @param min_junction_reads Integer. Minimum \code{total_junction_reads}.
#'   Default: \code{NULL}.
#' @param max_junction_reads Integer. Maximum \code{total_junction_reads}.
#'   Default: \code{NULL}.
#' @param min_pct_covered Numeric (0–100). Minimum \code{pct_events_covered}.
#'   Default: \code{NULL}.
#' @param custom_filters Named list of two-element numeric vectors
#'   \code{c(min, max)} applied to arbitrary \code{isoform_metadata} columns.
#'   Use \code{NA} for a one-sided bound. Default: \code{NULL}.
#' @param verbose Logical. Default: \code{TRUE}.
#'
#' @return The filtered \code{MatisseObject}.
#'
#' @seealso \code{\link{ComputeIsoformQC}}
#'
#' @rdname FilterCells
#' @export
setMethod("FilterCells", "MatisseObject",
          function(object,
                   min_junctions      = NULL,
                   max_junctions      = NULL,
                   min_junction_reads = NULL,
                   max_junction_reads = NULL,
                   min_pct_covered    = NULL,
                   custom_filters     = NULL,
                   verbose            = TRUE) {
  meta  <- object@isoform_metadata
  cells <- .get_cells(object)
  keep  <- rep(TRUE, length(cells))
  names(keep) <- cells

  apply_bound <- function(col, lo, hi) {
    if (!col %in% colnames(meta)) {
      rlang::warn(paste0("Column '", col,
                         "' not found in isoform_metadata. Run ComputeIsoformQC() first."))
      return()
    }
    vals <- meta[[col]]
    if (!is.null(lo) && !is.na(lo)) keep <<- keep & !is.na(vals) & (vals >= lo)
    if (!is.null(hi) && !is.na(hi)) keep <<- keep & !is.na(vals) & (vals <= hi)
  }

  apply_bound("n_junctions_detected",  min_junctions,      max_junctions)
  apply_bound("total_junction_reads",  min_junction_reads, max_junction_reads)
  apply_bound("pct_events_covered",    min_pct_covered,    NULL)

  if (!is.null(custom_filters)) {
    for (col in names(custom_filters)) {
      bounds <- custom_filters[[col]]
      apply_bound(col, bounds[1], bounds[2])
    }
  }

  n_removed <- sum(!keep)
  if (verbose) {
    cli::cli_alert_info(
      "Removing {n_removed} cells ({round(100*n_removed/length(keep),1)}%); ",
      "{sum(keep)} cells remain.")
  }

  object[cells[keep], ]
})

# ---------------------------------------------------------------------------
# FilterEvents
# ---------------------------------------------------------------------------

#' Filter splice events by coverage or variability
#'
#' Removes events that do not pass minimum coverage or variance thresholds.
#'
#' @param object A \code{MatisseObject} with a \code{"psi"} assay.
#' @param min_cells_covered Integer. Minimum number of cells in which the
#'   event must have a non-\code{NA} PSI value. Default: \code{10}.
#' @param min_psi_variance Numeric. Minimum variance of PSI across covered
#'   cells. Default: \code{NULL} (no variance filter).
#' @param verbose Logical. Default: \code{TRUE}.
#'
#' @return The filtered \code{MatisseObject}.
#'
#' @rdname FilterEvents
#' @export
setMethod("FilterEvents", "MatisseObject",
          function(object,
                   min_cells_covered = 10L,
                   min_psi_variance  = NULL,
                   verbose           = TRUE) {
  psi_cx <- GetPSI(object)  # cells x events
  if (is.null(psi_cx)) {
    rlang::abort("PSI matrix is NULL. Run CalculatePSI() first.")
  }

  n_covered <- .n_covered_per_event(psi_cx)
  keep      <- n_covered >= min_cells_covered

  if (!is.null(min_psi_variance)) {
    psi       <- .psi_to_dense_na(psi_cx)
    variances <- apply(psi, 2, stats::var, na.rm = TRUE)
    keep      <- keep & (variances >= min_psi_variance)
  }

  n_removed <- sum(!keep)
  if (verbose) {
    cli::cli_alert_info(
      "Removing {n_removed} events ({round(100*n_removed/length(keep),1)}%); ",
      "{sum(keep)} events remain.")
  }

  event_ids <- colnames(psi_cx)[keep]
  object[, event_ids]
})
