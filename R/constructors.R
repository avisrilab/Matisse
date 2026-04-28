#' @include MatisseObject-class.R
#' @include generics.R
NULL

#' Create a MatisseObject
#'
#' The single constructor for \code{\linkS4class{MatisseObject}}. Combines a
#' \code{Seurat} object with isoform-resolved splicing data and \emph{also}
#' computes PSI in one step (unless \code{defer_psi = TRUE}). The operating
#' mode is detected automatically from the inputs you supply:
#'
#' \itemize{
#'   \item \strong{Junction mode} (short-read): pass \code{junction_counts}
#'     plus splice events via \code{events}. Junction counts are stored as
#'     \code{Assay5("isoform")} inside the Seurat object; PSI values land in
#'     \code{Assay5("psi")}.
#'   \item \strong{Transcript mode} (long-read): pass \code{transcript_counts}
#'     plus splice events via \code{events} (typically a SUPPA2 IOE file).
#'     Transcript counts are stored as \code{Assay5("isoform")}; PSI values
#'     land in \code{Assay5("psi")}.
#' }
#'
#' Matisse intentionally does \strong{not} ship a heuristic event caller —
#' bring your own splice-event annotation (SUPPA2 IOE, rMATS table, or
#' equivalent).
#'
#' @param seurat A \code{Seurat} object. Required.
#' @param junction_counts A sparse matrix (dgCMatrix, cells x junctions) of
#'   raw per-junction read counts. Row names must match \code{colnames(seurat)}.
#'   Junction IDs (column names) ideally encode genomic coordinates
#'   (\code{chr-start-end-strand}, \code{chr:start-end:strand}, or
#'   \code{chr_start_end_strand}); Matisse parses these for sashimi plots.
#'   Triggers junction mode. Default: \code{NULL}.
#' @param transcript_counts A matrix or sparse matrix (transcripts x cells) of
#'   raw transcript-level counts. Stored as \code{Assay5("isoform")} in the
#'   Seurat object. Column names must overlap with \code{colnames(seurat)}.
#'   Triggers transcript mode. Default: \code{NULL}.
#' @param events Splice-event annotation. Either:
#'   \itemize{
#'     \item a character vector of paths to SUPPA2 \code{.ioe} file(s)
#'       (parsed internally), or
#'     \item a pre-built \code{data.frame} with columns \code{event_id},
#'       \code{gene_id}, \code{chr}, \code{strand}, \code{event_type},
#'       \code{inclusion_features}, \code{exclusion_features}. The
#'       \code{*_features} columns hold the IDs that support each role —
#'       junction IDs in junction mode, transcript IDs in transcript mode.
#'   }
#'   Required.
#' @param min_coverage Integer. Minimum total reads per cell per event for the
#'   PSI calculation step. Default: \code{5L}. Forwarded to
#'   \code{\link{CalculatePSI}} when \code{defer_psi = FALSE}.
#' @param defer_psi Logical. Skip the PSI calculation step at construction.
#'   The returned object will have raw counts only (no \code{"psi"} assay).
#'   Call \code{\link{CalculatePSI}} later. Default: \code{FALSE}.
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
#' # Junction mode (short-read): events as a data.frame
#' jxn <- make_junction_counts()
#' obj <- CreateMatisseObject(seurat = seu, junction_counts = jxn,
#'                            events = my_event_df)
#'
#' # Transcript mode (long-read): events as SUPPA2 IOE file paths
#' tx  <- make_transcript_counts()
#' obj <- CreateMatisseObject(seurat = seu, transcript_counts = tx,
#'                            events = "path/to/events.ioe")
#' }
#'
#' @export
CreateMatisseObject <- function(
    seurat,
    junction_counts   = NULL,
    transcript_counts = NULL,
    events            = NULL,
    min_coverage      = 5L,
    defer_psi         = FALSE,
    verbose           = TRUE
) {
  if (!inherits(seurat, "Seurat")) {
    rlang::abort("`seurat` must be a Seurat object.")
  }
  cells <- colnames(seurat)

  # Detect mode ------------------------------------------------------------
  has_junctions   <- !is.null(junction_counts)
  has_transcripts <- !is.null(transcript_counts)

  if (has_junctions && has_transcripts) {
    rlang::abort(paste0(
      "Provide either `junction_counts` (junction mode) or ",
      "`transcript_counts` (transcript mode), not both."))
  }

  input_mode <- if (has_transcripts) "transcript" else "junction"

  # --- resolve the events argument: paths -> parsed data.frame ------------
  events_path <- NA_character_
  event_data  <- NULL
  if (!is.null(events)) {
    if (is.character(events)) {
      missing_files <- events[!file.exists(events)]
      if (length(missing_files) > 0L) {
        rlang::abort(paste0("Event file(s) not found: ",
                            paste(missing_files, collapse = ", ")))
      }
      if (verbose) cli::cli_alert_info("Parsing {length(events)} event file(s)...")
      ioe         <- .parse_ioe_files(events)
      event_data  <- data.frame(
        event_id           = ioe$event_id,
        gene_id            = ioe$gene_id,
        chr                = ioe$chr,
        strand             = ioe$strand,
        event_type         = ioe$event_type,
        inclusion_features = ioe$inclusion_features,
        exclusion_features = ioe$exclusion_features,
        stringsAsFactors   = FALSE
      )
      events_path <- if (length(events) == 1L) {
        normalizePath(events, mustWork = FALSE)
      } else {
        vapply(events, normalizePath, character(1), mustWork = FALSE)
      }
    } else if (is.data.frame(events)) {
      required <- c("event_id", "gene_id", "chr", "strand", "event_type",
                    "inclusion_features", "exclusion_features")
      .check_required_columns(events, required, "events")
      event_data <- as.data.frame(events)
    } else {
      rlang::abort(paste0(
        "`events` must be a character vector of file paths or a data.frame; ",
        "got '", class(events)[1L], "'."))
    }
  }

  if (verbose) {
    cli::cli_alert_info(paste0(
      "Creating MatisseObject ({input_mode} mode): {length(cells)} cells, ",
      "{nrow(seurat)} gene features."))
  }

  # --- junction mode: store as Assay5("isoform") ---------------------------
  if (has_junctions) {
    junction_counts <- .validate_cell_matrix(junction_counts, cells,
                                             "junction_counts")
    seurat <- .add_isoform_assay_junction(seurat, junction_counts, verbose)
    # Auto-derive per-junction genomic coordinates from junction IDs and
    # write them to the isoform assay's meta.features. PlotSashimi reads
    # these in junction mode.
    seurat <- .write_junction_meta_features(seurat, verbose = verbose)
  }

  # --- transcript mode: subset seurat to cells present in tx counts --------
  if (has_transcripts) {
    common_cells <- intersect(colnames(transcript_counts), cells)
    if (length(common_cells) == 0) {
      rlang::abort(
        "No cell barcodes overlap between `transcript_counts` and `seurat`.")
    }
    if (length(common_cells) < length(cells)) {
      rlang::warn(paste0(
        "`transcript_counts` covers ", length(common_cells), "/",
        length(cells), " Seurat cells. Subsetting Seurat object."))
      seurat <- seurat[, common_cells]
    }
    seurat <- .add_isoform_assay_transcript(seurat, transcript_counts,
                                            common_cells, verbose)
  }

  # --- write nCount_isoform / nFeature_isoform to meta.data ----------------
  if (has_junctions || has_transcripts) {
    seurat <- .write_isoform_qc(seurat)
  }

  # --- stage event_data in @misc; CalculatePSI migrates to assay -----------
  obj_misc <- list()
  obj_misc[["event_data"]] <-
    if (!is.null(event_data) && nrow(event_data) > 0L) event_data
    else data.frame()
  obj_misc[["event_data_path"]] <- events_path

  version_str <- tryCatch(
    as.character(utils::packageVersion("Matisse")),
    error = function(e) "development"
  )

  obj <- methods::new(
    "MatisseObject",
    seurat     = seurat,
    input.mode = input_mode,
    version    = version_str,
    misc       = obj_misc
  )

  if (verbose) cli::cli_alert_success("MatisseObject created successfully.")

  methods::validObject(obj)

  # P5: fold PSI calculation into construction. Single-step UX — by the time
  # the constructor returns, the PSI assay exists and event annotation lives
  # in seurat[["psi"]][[]] (per P1). Skip when defer_psi=TRUE or when no
  # events were supplied (an isoform-only object with no PSI is still useful
  # for raw-count exploration).
  has_events <- !is.null(obj_misc[["event_data"]]) &&
                nrow(obj_misc[["event_data"]]) > 0L
  if (!defer_psi && has_events && (has_junctions || has_transcripts)) {
    obj <- CalculatePSI(obj, min_coverage = min_coverage, verbose = verbose)
  }

  obj
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Store junction counts as Assay5("isoform") in the Seurat object.
# Input: cells x junctions (Matisse convention); stored as junctions x cells.
.add_isoform_assay_junction <- function(seurat, jxn_counts, verbose) {
  if (ncol(jxn_counts) < 2L) {
    rlang::abort(paste0(
      "`junction_counts` must have at least 2 junctions ",
      "(Assay5 requires >=2 features)."))
  }
  # Transpose: cells x junctions -> junctions x cells for Seurat convention
  jxn_ec    <- Matrix::t(jxn_counts)
  iso_assay <- SeuratObject::CreateAssay5Object(counts = jxn_ec)
  seurat[["isoform"]] <- iso_assay
  if (verbose) {
    cli::cli_alert_info("Added 'isoform' assay: {ncol(jxn_counts)} junctions.")
  }
  seurat
}

# Add transcript counts as Assay5("isoform") in the Seurat object.
# Input: transcripts x cells (already in Seurat orientation).
.add_isoform_assay_transcript <- function(seurat, tx_counts, cells, verbose) {
  if (!is.matrix(tx_counts) && !inherits(tx_counts, "Matrix")) {
    rlang::abort(paste0(
      "`transcript_counts` must be a matrix or sparse Matrix ",
      "(transcripts x cells)."))
  }
  if (is.null(colnames(tx_counts))) {
    rlang::abort("`transcript_counts` must have column names (cell barcodes).")
  }
  if (is.null(rownames(tx_counts))) {
    rlang::abort("`transcript_counts` must have row names (transcript IDs).")
  }
  tx_sub    <- tx_counts[, cells, drop = FALSE]
  tx_csc    <- methods::as(tx_sub, "CsparseMatrix")
  iso_assay <- SeuratObject::CreateAssay5Object(counts = tx_csc)
  seurat[["isoform"]] <- iso_assay
  if (verbose) {
    cli::cli_alert_info(
      "Added 'isoform' assay: {nrow(tx_csc)} transcripts.")
  }
  seurat
}

# Auto-derive per-junction genomic coordinates from junction IDs and write
# to the "isoform" assay's feature metadata. Junction-mode users get sashimi
# coordinates "for free" without supplying junction_data; the parser tolerates
# multiple ID formats and produces NA coords (with a warning) for anything it
# can't parse — counts still work, only sashimi degrades.
.write_junction_meta_features <- function(seurat, verbose = TRUE) {
  iso <- seurat[["isoform"]]
  jxn_ids <- rownames(iso)
  parsed <- .parse_junction_names(jxn_ids)
  rownames(parsed) <- parsed$junction_id
  parsed$junction_id <- NULL  # rownames carry the ID
  iso[[]] <- parsed
  seurat[["isoform"]] <- iso
  if (verbose) {
    n_ok <- sum(!is.na(parsed$chr))
    cli::cli_alert_info(paste0(
      "Parsed coordinates for {n_ok}/{nrow(parsed)} junctions ",
      "into the 'isoform' assay's meta.features."))
  }
  seurat
}

# Write nCount_isoform / nFeature_isoform to seurat@meta.data from "isoform" assay.
.write_isoform_qc <- function(seurat) {
  iso_assay  <- seurat[["isoform"]]
  counts_ev  <- SeuratObject::GetAssayData(iso_assay, layer = "counts")
  meta_df <- data.frame(
    nCount_isoform   = as.integer(Matrix::colSums(counts_ev)),
    nFeature_isoform = as.integer(Matrix::colSums(counts_ev > 0)),
    row.names        = colnames(counts_ev)
  )
  SeuratObject::AddMetaData(seurat, meta_df)
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

  events <- data.frame(
    event_id           = event_ids,
    gene_id            = gene_ids,
    chr                = ioe$seqname,
    strand             = strand,
    event_type         = event_types,
    inclusion_features = vapply(inc_list, paste, character(1), collapse = ";"),
    exclusion_features = vapply(exc_list, paste, character(1), collapse = ";"),
    stringsAsFactors   = FALSE,
    row.names          = NULL
  )

  # The same event ID can appear in multiple IOE files (SUPPA2 generates one
  # file per event type, and the same splice site can be annotated in several).
  # Assay5 requires unique rownames, so deduplicate: keep the first occurrence.
  dup_mask <- duplicated(events$event_id)
  if (any(dup_mask)) {
    rlang::warn(paste0(
      sum(dup_mask), " duplicate event ID(s) removed after merging IOE files. ",
      "Keeping the first occurrence of each."
    ))
    events <- events[!dup_mask, , drop = FALSE]
  }

  events
}

# ---------------------------------------------------------------------------
# Internal: aggregate transcript counts to per-event inclusion/exclusion
# ---------------------------------------------------------------------------

.aggregate_transcript_counts <- function(tx_counts, events,
                                          min_coverage, cells) {
  tx_names  <- rownames(tx_counts)

  inc_lists <- strsplit(events$inclusion_features, ";", fixed = TRUE)
  exc_lists <- strsplit(events$exclusion_features, ";", fixed = TRUE)

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
