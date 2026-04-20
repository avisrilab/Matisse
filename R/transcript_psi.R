#' @include MatisseObject-class.R psi.R
NULL

# ---------------------------------------------------------------------------
# CreateMatisseObjectFromTranscripts
# ---------------------------------------------------------------------------

#' Create a MatisseObject from transcript-level counts
#'
#' Constructs a \code{MatisseObject} directly from a transcript-by-cell count
#' matrix (e.g. from Bagpiper or other transcript-aware quantifiers) by mapping
#' transcripts to splice events using SUPPA2 IOE annotation files.
#'
#' Transcript counts are stored as a \code{Assay5} named \code{"transcript"}
#' inside the Seurat object (ready for \code{\link{SCTransformTranscripts}}).
#' PSI values, inclusion counts, and exclusion counts are stored inside a
#' \code{ChromatinAssay} named \code{"psi"} with layers \code{"data"},
#' \code{"counts"}, and \code{"exclusion"}, respectively.
#'
#' For each splice event the PSI per cell is:
#' \deqn{PSI_{c,e} = \frac{\sum \text{inclusion transcript counts}}
#'                        {\sum \text{inclusion transcript counts} +
#'                         \sum \text{exclusion transcript counts}}}
#'
#' Cells where the total count (inclusion + exclusion) falls below
#' \code{min_coverage} are set to \code{NA}.
#'
#' @param seurat A \code{Seurat} object. Provides cell barcodes and gene-level
#'   metadata.
#' @param transcript_counts A matrix or sparse matrix (\strong{transcripts ×
#'   cells}) of transcript-level counts. Row names must be transcript IDs
#'   matching those in the IOE files. Column names must be cell barcodes.
#' @param ioe_files Character vector of paths to SUPPA2 \code{.ioe} files.
#'   Multiple files (one per event type: SE, SS, MX, RI, FL) can be supplied.
#' @param min_coverage Integer. Minimum total transcript counts per cell per
#'   event required to report a PSI value. Default: \code{5}.
#' @param verbose Logical. Print progress messages. Default: \code{TRUE}.
#'
#' @return A \code{\linkS4class{MatisseObject}} with \code{"transcript"} and
#'   \code{"psi"} assays populated inside the embedded Seurat object.
#'
#' @seealso \code{\link{CreateMatisseObject}}, \code{\link{CalculatePSI}},
#'   \code{\link{SCTransformTranscripts}}
#'
#' @export
CreateMatisseObjectFromTranscripts <- function(
    seurat,
    transcript_counts,
    ioe_files,
    min_coverage = 5L,
    verbose      = TRUE
) {
  if (!inherits(seurat, "Seurat")) {
    rlang::abort("`seurat` must be a Seurat object.")
  }
  cells <- colnames(seurat)

  if (!is.matrix(transcript_counts) && !inherits(transcript_counts, "Matrix")) {
    rlang::abort(paste0(
      "`transcript_counts` must be a matrix or sparse Matrix (transcripts x cells)."))
  }
  if (is.null(rownames(transcript_counts))) {
    rlang::abort("`transcript_counts` must have row names (transcript IDs).")
  }
  if (is.null(colnames(transcript_counts))) {
    rlang::abort("`transcript_counts` must have column names (cell barcodes).")
  }

  common_cells <- intersect(colnames(transcript_counts), cells)
  if (length(common_cells) == 0) {
    rlang::abort(paste0(
      "No cell barcodes overlap between `transcript_counts` and `seurat`."))
  }
  if (length(common_cells) < length(cells)) {
    rlang::warn(paste0(
      length(common_cells), "/", length(cells),
      " Seurat cells found in `transcript_counts`. Only overlapping cells kept."))
  }
  transcript_counts <- transcript_counts[, common_cells, drop = FALSE]

  if (missing(ioe_files) || length(ioe_files) == 0) {
    rlang::abort("`ioe_files` must be provided: paths to SUPPA2 .ioe files.")
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
      "Found {nrow(events)} events across \\
      {length(unique(events$event_type))} event type(s); \\
      mapping to {nrow(transcript_counts)} transcripts...")
  }

  result <- .aggregate_transcript_counts(
    tx_counts    = transcript_counts,
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

  version_str <- tryCatch(
    as.character(utils::packageVersion("Matisse")),
    error = function(e) "development"
  )

  seurat_sub <- seurat[, common_cells]

  # Store transcript counts as "transcript" Assay5 in the Seurat object
  tx_csc <- methods::as(transcript_counts, "CsparseMatrix")
  seurat_sub[["transcript"]] <- SeuratObject::CreateAssay5Object(counts = tx_csc)

  # Store PSI as "psi" ChromatinAssay (data/counts/exclusion layers)
  seurat_sub[["psi"]] <- .create_psi_chromatin_assay(
    psi_mat      = result$psi,
    inc_mat      = result$inclusion,
    exc_mat      = result$exclusion,
    event_data   = event_data,
    junction_data = NULL
  )

  obj <- methods::new(
    "MatisseObject",
    seurat           = seurat_sub,
    junction_counts  = NULL,
    event_data       = event_data,
    junction_data    = data.frame(),
    isoform_metadata = data.frame(row.names = common_cells),
    version          = version_str,
    misc             = list(source = "transcript_counts")
  )

  if (verbose) {
    pct <- round(100 * sum(.n_covered_per_event(result$psi)) /
                   (as.double(nrow(result$psi)) * ncol(result$psi)), 1)
    cli::cli_alert_success(
      "MatisseObject created from transcript counts: \\
      {length(common_cells)} cells, {nrow(events)} events, \\
      {pct}% PSI entries covered.")
  }

  methods::validObject(obj)
  obj
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
