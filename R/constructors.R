#' @include MatisseObject-class.R
#' @include generics.R
NULL

#' Create a MatisseObject
#'
#' The single constructor for \code{\linkS4class{MatisseObject}}. Combines a
#' \code{Seurat} object with isoform-resolved splicing data. The operating
#' mode is detected automatically from the inputs you supply:
#'
#' \itemize{
#'   \item \strong{Junction mode} (short-read): pass \code{junction_counts}.
#'     Junction counts are stored as \code{Assay5("junction")} inside the
#'     Seurat object. Call \code{\link{CalculatePSI}} afterwards to compute
#'     PSI values.
#'   \item \strong{Event mode} (long-read): pass \code{transcript_counts} and
#'     \code{ioe_files}. Transcript counts are stored as
#'     \code{Assay5("transcript")} and PSI is computed immediately from the
#'     supplied SUPPA2 \code{.ioe} event definitions, stored as
#'     \code{Assay5("psi")}.
#' }
#'
#' You can also pass \code{transcript_counts} alone (without \code{ioe_files})
#' to store the transcript assay without computing PSI -- for example, when you
#' want to run \code{\link[Seurat]{SCTransform}} on transcript-level counts and compute
#' PSI separately.
#'
#' @param seurat A \code{Seurat} object. Required.
#' @param junction_counts A sparse matrix (dgCMatrix, cells x junctions) of
#'   raw per-junction read counts. Row names must match \code{colnames(seurat)}.
#'   Triggers junction mode. Default: \code{NULL}.
#' @param transcript_counts A matrix or sparse matrix (transcripts x cells) of
#'   raw transcript-level counts. Stored as \code{Assay5("transcript")} in the
#'   Seurat object. Column names must overlap with \code{colnames(seurat)}.
#'   Default: \code{NULL}.
#' @param ioe_files Character vector of paths to SUPPA2 \code{.ioe} files.
#'   When supplied together with \code{transcript_counts}, triggers event mode
#'   and PSI is computed at construction. Default: \code{NULL}.
#' @param event_data A \code{data.frame} defining splice events (junction mode
#'   only). Required columns: \code{event_id}, \code{gene_id}, \code{chr},
#'   \code{strand}, \code{event_type}, \code{inclusion_junctions},
#'   \code{exclusion_junctions}. Default: \code{NULL}.
#' @param junction_data A \code{data.frame} of junction annotations. Required
#'   columns: \code{junction_id}, \code{chr}, \code{start}, \code{end},
#'   \code{strand}, \code{gene_id}. Default: \code{NULL}.
#' @param min_coverage Integer. Minimum total transcript counts per cell per
#'   event to report a PSI value (event mode only). Default: \code{5}.
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
#'
#' # Junction mode
#' jxn <- make_junction_counts()
#' obj <- CreateMatisseObject(seurat = seu, junction_counts = jxn)
#'
#' # Event mode
#' tx  <- make_transcript_counts()
#' obj <- CreateMatisseObject(seurat = seu, transcript_counts = tx,
#'                            ioe_files = "path/to/events.ioe")
#' }
#'
#' @export
CreateMatisseObject <- function(
    seurat,
    junction_counts   = NULL,
    transcript_counts = NULL,
    ioe_files         = NULL,
    event_data        = NULL,
    junction_data     = NULL,
    min_coverage      = 5L,
    verbose           = TRUE
) {
  if (!inherits(seurat, "Seurat")) {
    rlang::abort("`seurat` must be a Seurat object.")
  }
  cells <- colnames(seurat)

  # Detect mode ------------------------------------------------------------
  has_junctions   <- !is.null(junction_counts)
  has_ioe         <- !is.null(ioe_files) && length(ioe_files) > 0
  has_transcripts <- !is.null(transcript_counts)

  if (has_junctions && has_ioe) {
    rlang::abort(
      "Provide either `junction_counts` (junction mode) or `ioe_files` ",
      "(event mode), not both.")
  }

  mode <- if (has_ioe) "event" else "junction"

  if (verbose) {
    cli::cli_alert_info(
      "Creating MatisseObject ({mode} mode): {length(cells)} cells, ",
      "{nrow(seurat)} gene features.")
  }

  # --- validate and store junction_counts as Assay5("junction") ------------
  if (has_junctions) {
    junction_counts <- .validate_cell_matrix(junction_counts, cells,
                                              "junction_counts")
    seurat <- .add_junction_assay(seurat, junction_counts, cells, verbose)
  }

  # --- optional transcript counts -> "transcript" Assay5 -------------------
  if (has_transcripts) {
    seurat <- .add_transcript_assay(seurat, transcript_counts, cells, verbose)
  }

  # --- validate event_data -------------------------------------------------
  if (is.null(event_data)) {
    event_data <- data.frame()
  } else {
    required <- c("event_id", "gene_id", "chr", "strand",
                  "event_type", "inclusion_junctions", "exclusion_junctions")
    .check_required_columns(event_data, required, "event_data")
    event_data <- as.data.frame(event_data)
  }

  # --- validate junction_data ----------------------------------------------
  if (is.null(junction_data)) {
    junction_data <- data.frame()
  } else {
    required <- c("junction_id", "chr", "start", "end", "strand", "gene_id")
    .check_required_columns(junction_data, required, "junction_data")
    junction_data <- as.data.frame(junction_data)
  }

  # --- event mode: parse IOE files and compute PSI -------------------------
  if (has_ioe) {
    if (!has_transcripts) {
      rlang::abort(
        "`transcript_counts` must be supplied together with `ioe_files` ",
        "for event-mode construction.")
    }

    missing_files <- ioe_files[!file.exists(ioe_files)]
    if (length(missing_files) > 0) {
      rlang::abort(paste0("IOE file(s) not found: ",
                          paste(missing_files, collapse = ", ")))
    }

    if (verbose) cli::cli_alert_info("Parsing {length(ioe_files)} IOE file(s)...")
    events <- .parse_ioe_files(ioe_files)

    if (verbose) {
      cli::cli_alert_info(
        "Found {nrow(events)} events across ",
        "{length(unique(events$event_type))} event type(s); ",
        "mapping to {nrow(transcript_counts)} transcripts...")
    }

    common_cells <- intersect(colnames(transcript_counts), cells)
    result <- .aggregate_transcript_counts(
      tx_counts    = transcript_counts[, common_cells, drop = FALSE],
      events       = events,
      min_coverage = min_coverage,
      cells        = common_cells
    )

    event_data <- data.frame(
      event_id             = events$event_id,
      gene_id              = events$gene_id,
      chr                  = events$chr,
      strand               = events$strand,
      event_type           = events$event_type,
      inclusion_junctions  = events$inclusion_transcripts,
      exclusion_junctions  = events$exclusion_transcripts,
      stringsAsFactors     = FALSE
    )

    # Subset Seurat to cells that are in transcript_counts
    seurat <- seurat[, common_cells]

    psi_result <- .create_psi_assay(
      psi_mat = result$psi,
      inc_mat = result$inclusion,
      exc_mat = result$exclusion
    )
    seurat[["psi"]] <- psi_result$assay
    # Sync event_data$event_id with stored names (SeuratObject may sanitize)
    event_data$event_id <- psi_result$feature_names

    if (verbose) {
      pct <- round(100 * sum(.n_covered_per_event(result$psi)) /
                     (as.double(nrow(result$psi)) * ncol(result$psi)), 1)
      cli::cli_alert_success(
        "PSI computed from transcript counts: {nrow(events)} events, ",
        "{pct}% entries covered.")
    }
  }

  version_str <- tryCatch(
    as.character(utils::packageVersion("Matisse")),
    error = function(e) "development"
  )

  obj <- methods::new(
    "MatisseObject",
    seurat        = seurat,
    event_data    = event_data,
    junction_data = junction_data,
    mode          = mode,
    version       = version_str,
    misc          = list()
  )

  if (verbose) cli::cli_alert_success("MatisseObject created successfully.")

  methods::validObject(obj)
  obj
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Store junction counts as Assay5("junction") in the Seurat object.
# Input: cells x junctions (Matisse convention); stored as junctions x cells.
.add_junction_assay <- function(seurat, jxn_counts, cells, verbose) {
  if (ncol(jxn_counts) < 2L) {
    rlang::abort(
      "`junction_counts` must have at least 2 junctions ",
      "(Assay5 requires >=2 features).")
  }
  # Transpose: cells x junctions -> junctions x cells for Seurat convention
  jxn_ec  <- Matrix::t(jxn_counts)
  jxn_assay <- SeuratObject::CreateAssay5Object(counts = jxn_ec)
  seurat[["junction"]] <- jxn_assay
  if (verbose) {
    cli::cli_alert_info("Added 'junction' assay: {ncol(jxn_counts)} junctions.")
  }
  seurat
}

# Add a transcript count matrix as an Assay5 named "transcript" in the Seurat
# object. The matrix must be transcripts x cells.
.add_transcript_assay <- function(seurat, tx_counts, cells, verbose) {
  if (!is.matrix(tx_counts) && !inherits(tx_counts, "Matrix")) {
    rlang::abort(
      "`transcript_counts` must be a matrix or sparse Matrix (transcripts x cells).")
  }
  if (is.null(colnames(tx_counts))) {
    rlang::abort("`transcript_counts` must have column names (cell barcodes).")
  }
  if (is.null(rownames(tx_counts))) {
    rlang::abort("`transcript_counts` must have row names (transcript IDs).")
  }
  common <- intersect(colnames(tx_counts), cells)
  if (length(common) == 0) {
    rlang::abort(
      "No cell barcodes overlap between `transcript_counts` and `seurat`.")
  }
  if (length(common) < length(cells)) {
    rlang::warn(paste0(
      "`transcript_counts` covers ", length(common), "/", length(cells),
      " Seurat cells."))
  }
  tx_sub   <- tx_counts[, common, drop = FALSE]
  tx_csc   <- methods::as(tx_sub, "CsparseMatrix")
  tx_assay <- SeuratObject::CreateAssay5Object(counts = tx_csc)
  seurat[["transcript"]] <- tx_assay
  if (verbose) {
    cli::cli_alert_info(
      "Added 'transcript' assay: {nrow(tx_csc)} transcripts.")
  }
  seurat
}

# Build an Assay5 ("psi") from PSI, inclusion, and exclusion matrices.
# All inputs are cells x events (Matisse convention); stored as events x cells.
# Returns list($assay, $feature_names) -- stored names may differ from inputs
# if SeuratObject sanitized them (e.g. underscore -> dash).
.create_psi_assay <- function(psi_mat, inc_mat, exc_mat) {
  psi_ec <- Matrix::t(psi_mat)
  inc_ec <- Matrix::t(inc_mat)
  exc_ec <- Matrix::t(exc_mat)

  assay <- SeuratObject::CreateAssay5Object(counts = inc_ec)
  # SeuratObject may sanitize feature names; read back and align the other layers
  stored_names <- rownames(assay)
  rownames(psi_ec) <- stored_names
  rownames(exc_ec) <- stored_names

  assay <- SeuratObject::SetAssayData(assay, layer = "data",      new.data = psi_ec)
  assay <- SeuratObject::SetAssayData(assay, layer = "exclusion", new.data = exc_ec)
  list(assay = assay, feature_names = stored_names)
}

# Coerce to dgCMatrix and check row names align with cell barcodes.
.validate_cell_matrix <- function(mat, cells, name) {
  if (!inherits(mat, "Matrix")) {
    mat <- Matrix::Matrix(as.matrix(mat), sparse = TRUE)
  }
  if (is.null(rownames(mat))) {
    rlang::abort(paste0("`", name, "` must have row names matching cell barcodes."))
  }
  if (!identical(sort(rownames(mat)), sort(cells))) {
    common <- intersect(rownames(mat), cells)
    if (length(common) == 0) {
      rlang::abort(paste0(
        "No row names of `", name, "` match cell barcodes in `seurat`."))
    }
    if (length(common) < length(cells)) {
      rlang::warn(paste0(
        "`", name, "` covers ", length(common), "/", length(cells),
        " cells. Missing cells will have all-zero rows."))
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

# Abort if a data.frame is missing required column names.
.check_required_columns <- function(df, required, name) {
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    rlang::abort(paste0(
      "`", name, "` is missing required columns: ",
      paste(missing, collapse = ", ")))
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Internal: parse SUPPA2 IOE files
# ---------------------------------------------------------------------------

.parse_ioe_files <- function(ioe_files) {
  dfs <- lapply(ioe_files, function(f) {
    df <- utils::read.table(
      f, header = TRUE, sep = "\t",
      stringsAsFactors = FALSE, quote = "", comment.char = ""
    )
    if (ncol(df) < 4) {
      rlang::abort(paste0(
        "IOE file '", f, "' must have at least 4 columns."))
    }
    col2_has_semicolon <- any(grepl(";", df[[2L]], fixed = TRUE))
    if (!col2_has_semicolon && ncol(df) >= 5L &&
        any(grepl(";", df[[3L]], fixed = TRUE))) {
      df <- df[, c(1L, 3L, 4L, 5L)]
    }
    colnames(df)[1:4] <- c("seqname", "gene_event_id",
                           "inclusion_transcripts", "total_transcripts")
    df        <- df[, 1:4]
    df$.source <- f
    df
  })
  ioe <- do.call(rbind, dfs)

  parts     <- strsplit(ioe$gene_event_id, ";", fixed = TRUE)
  gene_ids  <- vapply(parts, function(x) x[1L], character(1))
  event_ids <- vapply(parts, function(x)
    if (length(x) >= 2L) x[2L] else NA_character_, character(1))

  if (any(is.na(event_ids))) {
    bad    <- which(is.na(event_ids))
    n_show <- min(5L, length(bad))
    examples <- paste(
      sprintf(
        "  row %d in '%s'\n    seqname: %s\n    gene_id: %s",
        bad[seq_len(n_show)],
        ioe$.source[bad[seq_len(n_show)]],
        ioe$seqname[bad[seq_len(n_show)]],
        ioe$gene_event_id[bad[seq_len(n_show)]]
      ),
      collapse = "\n"
    )
    rlang::abort(paste0(
      "Some rows in the IOE file(s) have a malformed gene_id column. ",
      "Expected format: 'gene_id;event_id'.\n",
      length(bad), " problematic row(s) found. First ", n_show, ":\n",
      examples))
  }

  ioe$.source <- NULL

  event_types <- sub(":.*", "", event_ids)
  strand      <- sub(".*:", "", event_ids)

  inc_list <- strsplit(ioe$inclusion_transcripts, ",", fixed = TRUE)
  tot_list <- strsplit(ioe$total_transcripts,     ",", fixed = TRUE)
  exc_list <- mapply(setdiff, tot_list, inc_list, SIMPLIFY = FALSE)

  data.frame(
    event_id              = event_ids,
    gene_id               = gene_ids,
    chr                   = ioe$seqname,
    strand                = strand,
    event_type            = event_types,
    inclusion_transcripts = vapply(inc_list, paste, character(1), collapse = ";"),
    exclusion_transcripts = vapply(exc_list, paste, character(1), collapse = ";"),
    stringsAsFactors      = FALSE,
    row.names             = NULL
  )
}

# ---------------------------------------------------------------------------
# Internal: aggregate transcript counts to per-event inclusion/exclusion
# ---------------------------------------------------------------------------

.aggregate_transcript_counts <- function(tx_counts, events,
                                          min_coverage, cells) {
  tx_names  <- rownames(tx_counts)

  inc_lists <- strsplit(events$inclusion_transcripts, ";", fixed = TRUE)
  exc_lists <- strsplit(events$exclusion_transcripts, ";", fixed = TRUE)

  A_inc <- .build_indicator_matrix(inc_lists, tx_names)
  A_exc <- .build_indicator_matrix(exc_lists, tx_names)
  colnames(A_inc) <- colnames(A_exc) <- events$event_id

  tx_t    <- Matrix::t(tx_counts)   # cells x transcripts
  inc_mat <- tx_t %*% A_inc          # cells x events
  exc_mat <- tx_t %*% A_exc
  dimnames(inc_mat) <- dimnames(exc_mat) <- list(cells, events$event_id)

  psi_mat <- .psi_from_sparse_counts(inc_mat, exc_mat, min_coverage)

  list(
    psi       = psi_mat,
    inclusion = Matrix::Matrix(round(inc_mat), sparse = TRUE),
    exclusion = Matrix::Matrix(round(exc_mat), sparse = TRUE)
  )
}
