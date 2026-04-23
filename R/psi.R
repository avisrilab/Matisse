#' @include MatisseObject-class.R
#' @include generics.R
NULL

# ---------------------------------------------------------------------------
# CalculatePSI — MatisseObject method
# ---------------------------------------------------------------------------

#' Calculate PSI matrix from junction counts
#'
#' Computes a Percent Spliced In (PSI) matrix for all splice events defined in
#' \code{event_data}. Only applies to objects in \strong{junction mode}; in
#' event mode PSI is computed at construction time.
#'
#' For each cell \eqn{c} and event \eqn{e}:
#'
#' \deqn{PSI_{c,e} = \frac{\sum \text{inclusion reads}}
#'                        {\sum \text{inclusion reads} +
#'                         \sum \text{exclusion reads}}}
#'
#' Results are stored inside the embedded Seurat object as a
#' \code{Assay5} named \code{"psi"}, with:
#' \itemize{
#'   \item \code{"data"} layer: PSI values in \eqn{[0,1]} (events × cells).
#'   \item \code{"counts"} layer: inclusion read counts (events × cells).
#'   \item \code{"exclusion"} layer: exclusion read counts (events × cells).
#' }
#' Entries where total coverage falls below \code{min_coverage} are set to
#' \code{NA} in the \code{"data"} layer.
#'
#' @param object A \code{\linkS4class{MatisseObject}} in junction mode, or a
#'   sparse matrix (cells × junctions).
#' @param events When \code{object} is a matrix: a \code{data.frame} with
#'   columns \code{event_id}, \code{inclusion_junctions}, and
#'   \code{exclusion_junctions}. When \code{object} is a
#'   \code{MatisseObject} this defaults to \code{GetEventData(object)}.
#' @param min_coverage Integer. Minimum total reads per cell per event to
#'   report a PSI value. Default: \code{5}.
#' @param na_fill Numeric. Replacement for low-coverage entries. Default:
#'   \code{NA_real_}.
#' @param verbose Logical. Print progress. Default: \code{TRUE}.
#'
#' @return
#' * \code{MatisseObject}: the input object with the \code{"psi"} assay
#'   populated inside the embedded Seurat object.
#' * matrix: a dense matrix (cells × events) of PSI values.
#'
#' @seealso \code{\link{ComputeIsoformQC}}, \code{\link{PlotHeatmap}}
#'
#' @rdname CalculatePSI
#' @export
setMethod("CalculatePSI", "MatisseObject",
          function(object, events = NULL, min_coverage = 5L,
                   na_fill = NA_real_, verbose = TRUE, ...) {
  # In event mode PSI is already computed at construction — warn and return
  if (object@mode == "event") {
    if (!is.null(.get_assay_safe(object@seurat, "psi"))) {
      rlang::warn(paste0(
        "This object is in event mode and already has PSI computed at ",
        "construction. CalculatePSI() only applies to junction-mode objects. ",
        "Returning the object unchanged."))
      return(object)
    }
  }

  jxn_counts <- GetJunctionCounts(object)
  if (is.null(jxn_counts)) {
    rlang::abort(paste0(
      "No junction assay found. ",
      "Provide junction_counts via CreateMatisseObject()."))
  }

  if (is.null(events)) events <- object@event_data
  if (nrow(events) == 0) {
    rlang::abort(
      "No splice events defined. Provide event_data either via ",
      "CreateMatisseObject() or the `events` argument.")
  }

  result <- .calculate_psi_matrix(
    jxn_counts   = jxn_counts,
    events       = events,
    min_coverage = min_coverage,
    na_fill      = na_fill,
    verbose      = verbose
  )

  # Store PSI, inclusion, and exclusion in an Assay5 named "psi"
  psi_result <- .create_psi_assay(
    psi_mat = result$psi,
    inc_mat = result$inclusion,
    exc_mat = result$exclusion
  )
  object@seurat[["psi"]] <- psi_result$assay
  # Sync event_data$event_id with the names actually stored
  if (nrow(object@event_data) > 0) {
    idx <- match(events$event_id, object@event_data$event_id)
    idx <- idx[!is.na(idx)]
    object@event_data$event_id[idx] <- psi_result$feature_names
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
                   na_fill = NA_real_, verbose = TRUE, ...) {
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
    na_fill      = na_fill,
    verbose      = verbose
  )
  .psi_to_dense_na(result$psi)
})

# ---------------------------------------------------------------------------
# Core PSI computation (internal)
# ---------------------------------------------------------------------------

.calculate_psi_matrix <- function(jxn_counts, events,
                                   min_coverage, na_fill, verbose) {
  .check_required_columns(
    events,
    c("event_id", "inclusion_junctions", "exclusion_junctions"),
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

  inc_lists <- lapply(strsplit(events$inclusion_junctions, ";", fixed = TRUE),
                      trimws)
  exc_lists <- lapply(strsplit(events$exclusion_junctions, ";", fixed = TRUE),
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

.parse_junction_list <- function(x) {
  if (is.na(x) || nchar(trimws(x)) == 0) return(character(0))
  trimws(strsplit(x, ";", fixed = TRUE)[[1]])
}

# ---------------------------------------------------------------------------
# SummarizePSI
# ---------------------------------------------------------------------------

#' Summarize PSI distribution across cells for each event
#'
#' Returns a summary table with per-event PSI statistics across all (or a
#' subset of) cells. Call this after \code{\link{CalculatePSI}} or after
#' creating an event-mode object with \code{\link{CreateMatisseObject}}.
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

  psi     <- .psi_to_dense_na(psi_cx)
  n_cov   <- .n_covered_per_event(psi_cx)
  no_data <- n_cov == 0L

  mean_psi   <- colMeans(psi, na.rm = TRUE)
  median_psi <- apply(psi, 2, stats::median, na.rm = TRUE)
  sd_psi     <- apply(psi, 2, stats::sd,     na.rm = TRUE)

  mean_psi[no_data]   <- NA_real_
  median_psi[no_data] <- NA_real_
  sd_psi[no_data]     <- NA_real_

  data.frame(
    event_id        = colnames(psi),
    mean_psi        = mean_psi,
    median_psi      = median_psi,
    sd_psi          = sd_psi,
    n_cells_covered = n_cov,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
}
