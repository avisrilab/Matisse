#' @include MatisseObject-class.R
#' @include generics.R
NULL

# ---------------------------------------------------------------------------
# CalculatePSI — MatisseObject method
# ---------------------------------------------------------------------------

#' Calculate PSI matrix from junction counts
#'
#' Computes a Percent Spliced In (PSI) matrix for all splice events defined in
#' \code{event_data}. For each cell \eqn{c} and event \eqn{e}:
#'
#' \deqn{PSI_{c,e} = \frac{\sum \text{inclusion reads}}
#'                        {\sum \text{inclusion reads} +
#'                         \sum \text{exclusion reads}}}
#'
#' Entries where the total coverage (inclusion + exclusion) falls below
#' \code{min_coverage} are set to \code{na_fill} (default \code{NA}).
#'
#' @param object A \code{\linkS4class{MatisseObject}} that has a non-\code{NULL}
#'   \code{junction_counts} slot, or a sparse matrix (cells x junctions).
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
#' * \code{MatisseObject}: the input object with \code{psi},
#'   \code{inclusion_counts}, and \code{exclusion_counts} slots populated.
#' * matrix: a dense matrix (cells x events) of PSI values.
#'
#' @seealso \code{\link{ComputeIsoformQC}}, \code{\link{PlotPSIHeatmap}}
#'
#' @examples
#' \dontrun{
#' # Build a small toy junction matrix
#' jxn_mat <- Matrix::sparseMatrix(
#'   i = c(1,1,2,2),
#'   j = c(1,2,2,3),
#'   x = c(10, 5, 8, 3),
#'   dims = c(3, 4),
#'   dimnames = list(
#'     paste0("Cell", 1:3),
#'     c("jxn1","jxn2","jxn3","jxn4")
#'   )
#' )
#' events <- data.frame(
#'   event_id             = "SE_gene1",
#'   gene_id              = "gene1",
#'   chr                  = "chr1",
#'   strand               = "+",
#'   event_type           = "SE",
#'   inclusion_junctions  = "jxn1;jxn2",
#'   exclusion_junctions  = "jxn3",
#'   stringsAsFactors     = FALSE
#' )
#' psi_mat <- CalculatePSI(jxn_mat, events, min_coverage = 3)
#' }
#'
#' @rdname CalculatePSI
#' @export
setMethod("CalculatePSI", "MatisseObject",
          function(object, events = NULL, min_coverage = 5L,
                   na_fill = NA_real_, verbose = TRUE) {
  if (is.null(object@junction_counts)) {
    if (!is.null(object@psi)) {
      rlang::warn(paste0(
        "This object was created from transcript counts and already has PSI computed. ",
        "CalculatePSI() only applies to junction-count objects (from CreateMatisseObject()). ",
        "Returning the object unchanged."))
      return(object)
    }
    rlang::abort(
      "junction_counts slot is NULL. Provide junction counts via CreateMatisseObject().")
  }
  if (is.null(events)) {
    events <- object@event_data
  }
  if (nrow(events) == 0) {
    rlang::abort(
      "No splice events defined. Provide event_data either via \\
       CreateMatisseObject() or the `events` argument.")
  }

  result <- .calculate_psi_matrix(
    jxn_counts    = object@junction_counts,
    events        = events,
    min_coverage  = min_coverage,
    na_fill       = na_fill,
    verbose       = verbose
  )

  object@psi              <- result$psi
  object@inclusion_counts <- result$inclusion
  object@exclusion_counts <- result$exclusion

  if (verbose) {
    pct <- round(100 * sum(.n_covered_per_event(result$psi)) /
                   (as.numeric(nrow(result$psi)) * ncol(result$psi)), 1)
    cli::cli_alert_success(
      "PSI calculated for {ncol(result$psi)} events; {pct}% entries covered.")
  }

  methods::validObject(object)
  object
})

# CalculatePSI for a raw matrix input
#' @rdname CalculatePSI
setMethod("CalculatePSI", "ANY",
          function(object, events, min_coverage = 5L,
                   na_fill = NA_real_, verbose = TRUE) {
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

#' Compute PSI matrix from a junction count matrix
#'
#' @param jxn_counts Sparse matrix (cells x junctions).
#' @param events data.frame with event definitions.
#' @param min_coverage Integer threshold.
#' @param na_fill Fill value for low-coverage entries.
#' @param verbose Logical.
#'
#' @return A list with elements \code{psi}, \code{inclusion}, \code{exclusion}.
#' @keywords internal
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

  # Parse all junction lists upfront (vectorised)
  inc_lists <- lapply(strsplit(events$inclusion_junctions, ";", fixed = TRUE),
                      trimws)
  exc_lists <- lapply(strsplit(events$exclusion_junctions, ";", fixed = TRUE),
                      trimws)

  # Build sparse indicator matrices: junctions x events
  A_inc <- .build_indicator_matrix(inc_lists, jxn_names)
  A_exc <- .build_indicator_matrix(exc_lists, jxn_names)
  colnames(A_inc) <- colnames(A_exc) <- events$event_id

  # Single matrix multiply replaces the per-event loop:
  #   (cells x junctions) %*% (junctions x events) → (cells x events)
  inc_mat <- jxn_counts %*% A_inc   # cells x events
  exc_mat <- jxn_counts %*% A_exc
  dimnames(inc_mat) <- dimnames(exc_mat) <- list(cells, events$event_id)

  psi_mat <- .psi_from_sparse_counts(inc_mat, exc_mat, min_coverage)

  list(
    psi       = psi_mat,
    inclusion = Matrix::Matrix(round(inc_mat), sparse = TRUE),
    exclusion = Matrix::Matrix(round(exc_mat), sparse = TRUE)
  )
}

# Split a semicolon-delimited junction string into a character vector
.parse_junction_list <- function(x) {
  if (is.na(x) || nchar(trimws(x)) == 0) return(character(0))
  trimws(strsplit(x, ";", fixed = TRUE)[[1]])
}

# ---------------------------------------------------------------------------
# Convenience: get per-event PSI summary statistics
# ---------------------------------------------------------------------------

#' Summarize PSI distribution across cells for each event
#'
#' @param object A \code{MatisseObject} with a non-\code{NULL} \code{psi} slot.
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
  if (is.null(object@psi)) {
    rlang::abort("PSI matrix is empty. Run CalculatePSI() first.")
  }

  psi_sp <- object@psi
  if (!is.null(cells)) {
    psi_sp <- psi_sp[cells, , drop = FALSE]
  }
  # Convert to dense NA matrix for na.rm-aware statistics.
  # Absent entries in the sparse matrix = not covered = NA.
  psi <- .psi_to_dense_na(psi_sp)

  data.frame(
    event_id        = colnames(psi),
    mean_psi        = colMeans(psi, na.rm = TRUE),
    median_psi      = apply(psi, 2, stats::median, na.rm = TRUE),
    sd_psi          = apply(psi, 2, stats::sd,     na.rm = TRUE),
    n_cells_covered = .n_covered_per_event(psi_sp),
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
}
