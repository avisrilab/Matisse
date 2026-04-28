#' @include MatisseObject-class.R
#' @include generics.R
NULL

# ---------------------------------------------------------------------------
# FilterCells
# ---------------------------------------------------------------------------

#' Filter cells by QC thresholds
#'
#' Removes cells that do not pass the specified thresholds on QC columns in
#' \code{MatisseMeta(object)} (the Seurat \code{meta.data}). The columns
#' \code{nCount_isoform} and \code{nFeature_isoform} are written at construction
#' by \code{\link{CreateMatisseObject}}; \code{nPercent_isoform} is written by
#' \code{\link{CalculatePSI}}.
#'
#' @param object A \code{MatisseObject}.
#' @param min_features_isoform Integer. Minimum \code{nFeature_isoform} (number
#'   of junctions or transcripts with at least one read). Default: \code{NULL}.
#' @param max_features_isoform Integer. Maximum \code{nFeature_isoform}.
#'   Default: \code{NULL}.
#' @param min_counts_isoform Integer. Minimum \code{nCount_isoform} (total reads
#'   in the isoform assay). Default: \code{NULL}.
#' @param max_counts_isoform Integer. Maximum \code{nCount_isoform}.
#'   Default: \code{NULL}.
#' @param min_pct_isoform Numeric (0-100). Minimum \code{nPercent_isoform}
#'   (percentage of splice events with a non-NA PSI value). Requires
#'   \code{\link{CalculatePSI}} to have been run. Default: \code{NULL}.
#' @param custom_filters Named list of two-element numeric vectors
#'   \code{c(min, max)} applied to arbitrary metadata columns.
#'   Use \code{NA} for a one-sided bound. Default: \code{NULL}.
#' @param verbose Logical. Default: \code{TRUE}.
#' @param ... Ignored; included for S4 generic compatibility.
#'
#' @return The filtered \code{MatisseObject}.
#'
#' @seealso \code{\link{FilterEvents}}, \code{\link{CalculatePSI}}
#'
#' @rdname FilterCells
#' @export
setMethod("FilterCells", "MatisseObject",
          function(object,
                   min_features_isoform = NULL,
                   max_features_isoform = NULL,
                   min_counts_isoform   = NULL,
                   max_counts_isoform   = NULL,
                   min_pct_isoform      = NULL,
                   custom_filters       = NULL,
                   verbose              = TRUE, ...) {
  meta  <- MatisseMeta(object)
  cells <- .get_cells(object)
  keep  <- rep(TRUE, length(cells))
  names(keep) <- cells

  apply_bound <- function(col, lo, hi) {
    if (!col %in% colnames(meta)) {
      rlang::warn(paste0("Column '", col, "' not found in metadata. ",
                         "Run CreateMatisseObject() or CalculatePSI() first."))
      return()
    }
    vals <- meta[[col]]
    if (!is.null(lo) && !is.na(lo)) keep <<- keep & !is.na(vals) & (vals >= lo)
    if (!is.null(hi) && !is.na(hi)) keep <<- keep & !is.na(vals) & (vals <= hi)
  }

  if (!is.null(min_features_isoform) || !is.null(max_features_isoform))
    apply_bound("nFeature_isoform", min_features_isoform, max_features_isoform)
  if (!is.null(min_counts_isoform) || !is.null(max_counts_isoform))
    apply_bound("nCount_isoform", min_counts_isoform, max_counts_isoform)
  if (!is.null(min_pct_isoform))
    apply_bound("nPercent_isoform", min_pct_isoform, NULL)

  if (!is.null(custom_filters)) {
    for (col in names(custom_filters)) {
      bounds <- custom_filters[[col]]
      apply_bound(col, bounds[1], bounds[2])
    }
  }

  n_removed <- sum(!keep)
  if (verbose) {
    cli::cli_alert_info(paste0(
      "Removing {n_removed} cells ({round(100*n_removed/length(keep),1)}%); ",
      "{sum(keep)} cells remain."))
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
#' @param ... Ignored; included for S4 generic compatibility.
#'
#' @return The filtered \code{MatisseObject}.
#'
#' @seealso \code{\link{FilterCells}}
#'
#' @rdname FilterEvents
#' @export
setMethod("FilterEvents", "MatisseObject",
          function(object,
                   min_cells_covered = 10L,
                   min_psi_variance  = NULL,
                   verbose           = TRUE, ...) {
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
    cli::cli_alert_info(paste0(
      "Removing {n_removed} events ({round(100*n_removed/length(keep),1)}%); ",
      "{sum(keep)} events remain."))
  }

  event_ids <- colnames(psi_cx)[keep]
  object[, event_ids]
})
