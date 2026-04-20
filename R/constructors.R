#' @include MatisseObject-class.R
#' @include generics.R
NULL

#' Create a MatisseObject
#'
#' The primary constructor for \code{\linkS4class{MatisseObject}}. Combines a
#' \code{Seurat} object with optional isoform-resolved data layers. When
#' \code{transcript_counts} is supplied the matrix is stored as a \code{Assay5}
#' named \code{"transcript"} inside the embedded Seurat object, ready for
#' downstream SCTransform normalisation.
#'
#' @param seurat A \code{Seurat} object. Required.
#' @param junction_counts A sparse matrix (dgCMatrix, cells × junctions) of
#'   raw per-junction read counts. Row names must match
#'   \code{colnames(seurat)}. Default: \code{NULL}.
#' @param transcript_counts A matrix or sparse matrix (transcripts × cells) of
#'   raw transcript-level counts. When supplied it is added to the Seurat
#'   object as a \code{Assay5} named \code{"transcript"}. Column names must
#'   overlap with \code{colnames(seurat)}. Default: \code{NULL}.
#' @param event_data A \code{data.frame} defining splice events. Required
#'   columns: \code{event_id}, \code{gene_id}, \code{chr}, \code{strand},
#'   \code{event_type}, \code{inclusion_junctions},
#'   \code{exclusion_junctions}. Default: \code{NULL}.
#' @param junction_data A \code{data.frame} of junction annotations. Required
#'   columns: \code{junction_id}, \code{chr}, \code{start}, \code{end},
#'   \code{strand}, \code{gene_id}. Default: \code{NULL}.
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
#' obj <- CreateMatisseObject(seurat = seu)
#' obj
#' }
#'
#' @export
CreateMatisseObject <- function(
    seurat,
    junction_counts   = NULL,
    transcript_counts = NULL,
    event_data        = NULL,
    junction_data     = NULL,
    verbose           = TRUE
) {
  if (!inherits(seurat, "Seurat")) {
    rlang::abort("`seurat` must be a Seurat object.")
  }
  cells <- colnames(seurat)

  if (verbose) {
    cli::cli_alert_info(
      "Creating MatisseObject: {length(cells)} cells, \\
       {nrow(seurat)} gene features.")
  }

  # --- validate junction_counts -----------------------------------------------
  if (!is.null(junction_counts)) {
    junction_counts <- .validate_cell_matrix(junction_counts, cells,
                                             "junction_counts")
  }

  # --- optional transcript counts → "transcript" Assay5 ----------------------
  if (!is.null(transcript_counts)) {
    seurat <- .add_transcript_assay(seurat, transcript_counts, cells, verbose)
  }

  # --- validate event_data ----------------------------------------------------
  if (is.null(event_data)) {
    event_data <- data.frame()
  } else {
    required <- c("event_id", "gene_id", "chr", "strand",
                  "event_type", "inclusion_junctions", "exclusion_junctions")
    .check_required_columns(event_data, required, "event_data")
    event_data <- as.data.frame(event_data)
  }

  # --- validate junction_data -------------------------------------------------
  if (is.null(junction_data)) {
    junction_data <- data.frame()
  } else {
    required <- c("junction_id", "chr", "start", "end", "strand", "gene_id")
    .check_required_columns(junction_data, required, "junction_data")
    junction_data <- as.data.frame(junction_data)
  }

  isoform_metadata <- data.frame(row.names = cells)

  version_str <- tryCatch(
    as.character(utils::packageVersion("Matisse")),
    error = function(e) "development"
  )

  obj <- methods::new(
    "MatisseObject",
    seurat           = seurat,
    junction_counts  = junction_counts,
    event_data       = event_data,
    junction_data    = junction_data,
    isoform_metadata = isoform_metadata,
    version          = version_str,
    misc             = list()
  )

  if (verbose) cli::cli_alert_success("MatisseObject created successfully.")

  methods::validObject(obj)
  obj
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Add a transcript count matrix as an Assay5 named "transcript" in the Seurat
# object. The matrix must be transcripts x cells.
.add_transcript_assay <- function(seurat, tx_counts, cells, verbose) {
  if (!is.matrix(tx_counts) && !inherits(tx_counts, "Matrix")) {
    rlang::abort("`transcript_counts` must be a matrix or sparse Matrix (transcripts x cells).")
  }
  if (is.null(colnames(tx_counts))) {
    rlang::abort("`transcript_counts` must have column names (cell barcodes).")
  }
  if (is.null(rownames(tx_counts))) {
    rlang::abort("`transcript_counts` must have row names (transcript IDs).")
  }
  common <- intersect(colnames(tx_counts), cells)
  if (length(common) == 0) {
    rlang::abort("No cell barcodes overlap between `transcript_counts` and `seurat`.")
  }
  if (length(common) < length(cells)) {
    rlang::warn(paste0(
      "`transcript_counts` covers ", length(common), "/", length(cells),
      " Seurat cells."))
  }
  tx_sub  <- tx_counts[, common, drop = FALSE]
  tx_csc  <- methods::as(tx_sub, "CsparseMatrix")
  tx_assay <- SeuratObject::CreateAssay5Object(counts = tx_csc)
  seurat[["transcript"]] <- tx_assay
  if (verbose) {
    cli::cli_alert_info(
      "Added 'transcript' assay: {nrow(tx_csc)} transcripts.")
  }
  seurat
}

# Build a ChromatinAssay ("psi") from PSI, inclusion, and exclusion matrices.
# All matrices are in cells x events (Matisse) convention; we transpose to
# events x cells (Seurat) before creating the assay.
.create_psi_chromatin_assay <- function(psi_mat, inc_mat, exc_mat,
                                         event_data, junction_data = NULL) {
  # Transpose: cells x events → events x cells (Seurat convention)
  psi_ec <- Matrix::t(psi_mat)
  inc_ec <- Matrix::t(inc_mat)
  exc_ec <- Matrix::t(exc_mat)

  gr <- .make_event_granges(event_data, junction_data)

  assay <- Signac::CreateChromatinAssay(
    counts  = inc_ec,
    data    = psi_ec,
    ranges  = gr
  )
  # Store exclusion counts as an additional named layer (Seurat v5 / Assay5)
  assay <- SeuratObject::SetAssayData(assay, layer = "exclusion",
                                       new.data = exc_ec)
  assay
}

# Build a GRanges object for splice events, using junction_data coordinates
# where available and placeholder ranges (start=1, width=100) otherwise.
.make_event_granges <- function(event_data, junction_data = NULL) {
  n        <- nrow(event_data)
  seqnames <- event_data$chr
  strand   <- event_data$strand
  starts   <- rep(1L, n)
  ends     <- rep(100L, n)

  if (!is.null(junction_data) && nrow(junction_data) > 0 &&
      all(c("junction_id", "start", "end") %in% colnames(junction_data))) {
    jxn_s <- stats::setNames(junction_data$start, junction_data$junction_id)
    jxn_e <- stats::setNames(junction_data$end,   junction_data$junction_id)

    for (i in seq_len(n)) {
      jxn_ids <- unlist(strsplit(
        paste(event_data$inclusion_junctions[i],
              event_data$exclusion_junctions[i], sep = ";"),
        ";", fixed = TRUE
      ))
      jxn_ids <- jxn_ids[nchar(trimws(jxn_ids)) > 0]
      s <- stats::na.omit(jxn_s[jxn_ids])
      e <- stats::na.omit(jxn_e[jxn_ids])
      if (length(s) > 0) {
        starts[i] <- min(as.integer(s))
        ends[i]   <- max(as.integer(e))
      }
    }
  }

  gr <- GenomicRanges::GRanges(
    seqnames = seqnames,
    ranges   = IRanges::IRanges(start = starts, end = ends),
    strand   = strand
  )
  names(gr) <- event_data$event_id
  gr
}

# Coerce to dgCMatrix and check row names align with cell barcodes
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

# Abort if a data.frame is missing required column names
.check_required_columns <- function(df, required, name) {
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    rlang::abort(paste0(
      "`", name, "` is missing required columns: ",
      paste(missing, collapse = ", ")))
  }
  invisible(TRUE)
}
