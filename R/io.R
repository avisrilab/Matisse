#' @include MatisseObject-class.R
NULL

# ---------------------------------------------------------------------------
# STARsolo splice-junction input
# ---------------------------------------------------------------------------

#' Read a STARsolo splice-junction (SJ) matrix
#'
#' Loads the per-cell splice-junction count matrix produced by STARsolo with
#' \code{--soloFeatures SJ} and relabels its junctions into the
#' \code{chr-start-end-strand} ID form Matisse expects for junction mode
#' (the first pattern recognised by the internal junction-ID parser, so the
#' counts are ready for \code{\link{CreateMatisseObject}} and sashimi plots).
#'
#' STARsolo writes the SJ matrix as \emph{junctions x cells} with a 9-column
#' \code{features.tsv} (\code{chr}, intron start, intron end, strand code,
#' splice motif, ...). The intron start/end are the first and last intronic
#' base; the strand code is \code{0} (undefined), \code{1} (\code{+}) or
#' \code{2} (\code{-}). Only a \code{raw} matrix is emitted for SJ, so pass
#' the filtered cell barcodes (e.g. from \code{Gene/filtered}) via
#' \code{cells} to subset.
#'
#' @param sj_dir Directory holding \code{matrix.mtx}, \code{barcodes.tsv} and
#'   \code{features.tsv} (optionally \code{.gz}). Typically a STARsolo
#'   \code{Solo.out/SJ/raw} directory.
#' @param cells Optional character vector of cell barcodes to keep (e.g.
#'   \code{colnames(seurat)}). Barcodes absent from the SJ matrix are dropped
#'   with a warning. Default: \code{NULL} (keep all).
#' @param strand_map Named character vector mapping STAR strand codes to
#'   symbols. Default: \code{c("0" = "*", "1" = "+", "2" = "-")}.
#' @param infer_strand Logical. When the STAR strand code is \code{0}
#'   (undefined) but the intron-motif column is canonical, recover the strand
#'   from the motif (\code{1/3/5 -> +}, \code{2/4/6 -> -}; \code{0} stays
#'   undefined). Roughly one in six STARsolo junctions are strand-undefined;
#'   without this they can never match the \code{+}/\code{-} junction IDs in a
#'   SUPPA2 event table. Default: \code{TRUE}.
#' @param verbose Logical. Print progress. Default: \code{TRUE}.
#'
#' @return A \code{dgCMatrix} of raw junction counts, \emph{cells x junctions}
#'   (Matisse convention), with junction IDs as column names. Feed directly to
#'   \code{CreateMatisseObject(junction_counts = ...)}.
#'
#' @seealso \code{\link{BuildJunctionEvents}} to turn SUPPA2 events into a
#'   junction-ID event table that matches these column names.
#'
#' @examples
#' \dontrun{
#' jxn <- ReadSTARsoloSJ("Solo.out/SJ/raw", cells = colnames(seu))
#' obj <- CreateMatisseObject(seu, junction_counts = jxn,
#'                            events = BuildJunctionEvents(ioe_files))
#' }
#'
#' @export
ReadSTARsoloSJ <- function(sj_dir, cells = NULL,
                           strand_map = c("0" = "*", "1" = "+", "2" = "-"),
                           infer_strand = TRUE,
                           verbose = TRUE) {
  if (!dir.exists(sj_dir)) {
    rlang::abort(paste0("`sj_dir` does not exist: ", sj_dir))
  }
  mtx_path  <- .resolve_starsolo_file(sj_dir, "matrix.mtx")
  bc_path   <- .resolve_starsolo_file(sj_dir, "barcodes.tsv")
  feat_path <- .resolve_starsolo_file(sj_dir, "features.tsv")

  if (verbose) cli::cli_alert_info("Reading STARsolo SJ features...")
  feat <- utils::read.delim(feat_path, header = FALSE, sep = "\t",
                            colClasses = "character", quote = "",
                            comment.char = "")
  if (ncol(feat) < 4L) {
    rlang::abort(paste0(
      "STARsolo SJ features.tsv must have >= 4 columns ",
      "(chr, intron start, intron end, strand code); got ", ncol(feat), "."))
  }
  strand_sym <- unname(strand_map[feat[[4L]]])
  strand_sym[is.na(strand_sym)] <- "*"
  if (infer_strand && ncol(feat) >= 5L) {
    undef  <- feat[[4L]] == "0"
    motif  <- feat[[5L]]
    plus_m  <- undef & motif %in% c("1", "3", "5")
    minus_m <- undef & motif %in% c("2", "4", "6")
    strand_sym[plus_m]  <- "+"
    strand_sym[minus_m] <- "-"
    if (verbose && (any(plus_m) || any(minus_m))) {
      cli::cli_alert_info(paste0(
        "Recovered strand from intron motif for ",
        "{sum(plus_m) + sum(minus_m)} strand-undefined junction(s)."))
    }
  }
  jxn_ids <- paste(feat[[1L]], feat[[2L]], feat[[3L]], strand_sym, sep = "-")

  barcodes <- readLines(bc_path)

  if (verbose) cli::cli_alert_info("Reading STARsolo SJ matrix...")
  m <- Matrix::readMM(mtx_path)            # junctions x cells
  if (nrow(m) != length(jxn_ids) || ncol(m) != length(barcodes)) {
    rlang::abort(paste0(
      "STARsolo SJ matrix dimensions (", nrow(m), " x ", ncol(m),
      ") do not match features (", length(jxn_ids), ") / barcodes (",
      length(barcodes), ")."))
  }
  rownames(m) <- jxn_ids
  colnames(m) <- barcodes

  if (!is.null(cells)) {
    keep <- intersect(cells, barcodes)
    if (length(keep) == 0L) {
      rlang::abort(paste0(
        "None of the requested `cells` are present in the SJ barcodes. ",
        "Example requested: ", paste(utils::head(cells, 3L), collapse = ", "),
        "; example SJ: ", paste(utils::head(barcodes, 3L), collapse = ", ")))
    }
    if (length(keep) < length(cells) && verbose) {
      cli::cli_alert_warning(paste0(
        "{length(keep)}/{length(cells)} requested cells found in the SJ ",
        "matrix; the rest had no junction reads and are dropped."))
    }
    m <- m[, keep, drop = FALSE]
  }

  dup <- duplicated(rownames(m))
  if (any(dup)) {
    rlang::warn(paste0(
      sum(dup), " duplicate junction ID(s) in STARsolo features (typically ",
      "strand-ambiguous calls collapsing onto the same coordinates); ",
      "keeping the first occurrence."))
    m <- m[!dup, , drop = FALSE]
  }

  out <- Matrix::t(m)                       # cells x junctions
  out <- methods::as(out, "CsparseMatrix")
  if (verbose) {
    cli::cli_alert_success(paste0(
      "STARsolo SJ: {nrow(out)} cells x {ncol(out)} junctions."))
  }
  out
}

# Resolve a STARsolo output file, tolerating an optional .gz suffix.
.resolve_starsolo_file <- function(dir, base) {
  cand <- file.path(dir, c(base, paste0(base, ".gz")))
  hit  <- cand[file.exists(cand)]
  if (length(hit) == 0L) {
    rlang::abort(paste0(
      "Could not find '", base, "' (or '", base, ".gz') in ", dir))
  }
  hit[1L]
}

# ---------------------------------------------------------------------------
# Long-read / transcript-quantifier input
# ---------------------------------------------------------------------------

#' Read a long-read transcript count matrix
#'
#' Loads a per-cell transcript count matrix from a 10x-style MatrixMarket
#' triplet (\code{matrix.mtx}, \code{barcodes.tsv}, \code{features.tsv},
#' optionally \code{.gz}) as produced by long-read / isoform quantifiers such
#' as \strong{Bagpiper}, \strong{FLAMES}, or \strong{LIQA}, and returns it in
#' the \emph{transcripts x cells} orientation with dimnames that
#' \code{\link{CreateMatisseObject}}'s transcript mode expects.
#'
#' Quantifiers disagree on matrix orientation (Bagpiper writes
#' \emph{cells x transcripts}; others write \emph{transcripts x cells}), so
#' \code{orientation = "auto"} infers it by matching the two matrix dimensions
#' against the feature and barcode counts and transposes as needed.
#'
#' \strong{Transcript-version landmine.} SUPPA2 \code{.ioe} events and the
#' quantifier may disagree on whether transcript IDs carry a version suffix
#' (\code{ENST...\.3}). PSI aggregation
#' (\code{.build_indicator_matrix}) matches IDs only after a \code{_}->\code{-}
#' sanitisation; it never strips versions, so a mismatch silently zeroes PSI.
#' \code{strip_version = TRUE} removes a trailing \code{.<digits>} from
#' transcript IDs (summing counts of any IDs that collapse) to defend against
#' this.
#'
#' @param mtx_dir Directory holding \code{matrix.mtx}, \code{barcodes.tsv} and
#'   \code{features.tsv} (optionally \code{.gz}).
#' @param cells Optional character vector of cell barcodes to keep (e.g.
#'   filtered barcodes). Barcodes absent from the matrix are dropped with a
#'   warning. Default: \code{NULL} (keep all).
#' @param strip_version Logical. Strip a trailing \code{.<digits>} version
#'   suffix from transcript IDs; counts of IDs that collapse to the same
#'   stripped ID are summed. Default: \code{TRUE}.
#' @param orientation One of \code{"auto"} (infer from dimensions; default),
#'   \code{"cells_x_tx"} (source is cells x transcripts; transpose), or
#'   \code{"tx_x_cells"} (source is already transcripts x cells).
#' @param verbose Logical. Print progress. Default: \code{TRUE}.
#'
#' @return A \code{dgCMatrix} of transcript counts, \emph{transcripts x cells},
#'   transcript IDs as row names and barcodes as column names. Feed directly
#'   to \code{CreateMatisseObject(transcript_counts = ...)}.
#'
#' @seealso \code{\link{ReadSTARsoloSJ}} for the short-read (junction) path;
#'   \code{\link{CreateMatisseObject}}.
#'
#' @examples
#' \dontrun{
#' txc <- ReadTranscriptMatrix("bagpiper/mats", cells = colnames(seu))
#' obj <- CreateMatisseObject(seu, transcript_counts = txc,
#'                            events = ioe_files)
#' }
#'
#' @export
ReadTranscriptMatrix <- function(mtx_dir, cells = NULL, strip_version = TRUE,
                                  orientation = c("auto", "cells_x_tx",
                                                  "tx_x_cells"),
                                  verbose = TRUE) {
  orientation <- match.arg(orientation)
  if (!dir.exists(mtx_dir)) {
    rlang::abort(paste0("`mtx_dir` does not exist: ", mtx_dir))
  }
  mtx_path  <- .resolve_starsolo_file(mtx_dir, "matrix.mtx")
  bc_path   <- .resolve_starsolo_file(mtx_dir, "barcodes.tsv")
  feat_path <- .resolve_starsolo_file(mtx_dir, "features.tsv")

  if (verbose) cli::cli_alert_info("Reading transcript features / barcodes...")
  # First whitespace/tab-delimited field; tolerates single-column files.
  tx_ids   <- sub("[\t ].*$", "", readLines(gzfile(feat_path, "rt")))
  barcodes <- sub("[\t ].*$", "", readLines(gzfile(bc_path, "rt")))
  nf <- length(tx_ids)
  nb <- length(barcodes)

  if (verbose) cli::cli_alert_info("Reading transcript matrix...")
  m  <- Matrix::readMM(gzfile(mtx_path))
  dr <- nrow(m)
  dc <- ncol(m)

  tx_rows   <- dr == nf && dc == nb
  cell_rows <- dr == nb && dc == nf
  ori <- switch(orientation,
    tx_x_cells = "tx_x_cells",
    cells_x_tx = "cells_x_tx",
    auto = if (tx_rows && !cell_rows) "tx_x_cells"
           else if (cell_rows && !tx_rows) "cells_x_tx"
           else if (tx_rows && cell_rows) {
             rlang::warn(paste0(
               "Matrix is square in features/barcodes; assuming ",
               "transcripts x cells. Set `orientation` explicitly if wrong."))
             "tx_x_cells"
           } else NA_character_)
  if (is.na(ori)) {
    rlang::abort(paste0(
      "Matrix dimensions (", dr, " x ", dc, ") match neither ",
      "transcripts x cells (", nf, " x ", nb, ") nor cells x transcripts (",
      nb, " x ", nf, "). Check the triplet files."))
  }

  if (ori == "cells_x_tx") {
    if (!(dr == nb && dc == nf)) {
      rlang::abort(paste0(
        "orientation='cells_x_tx' but matrix (", dr, " x ", dc,
        ") != barcodes x features (", nb, " x ", nf, ")."))
    }
    rownames(m) <- barcodes
    colnames(m) <- tx_ids
    m <- Matrix::t(m)
  } else {
    if (!(dr == nf && dc == nb)) {
      rlang::abort(paste0(
        "orientation='tx_x_cells' but matrix (", dr, " x ", dc,
        ") != features x barcodes (", nf, " x ", nb, ")."))
    }
    rownames(m) <- tx_ids
    colnames(m) <- barcodes
  }
  # m is now transcripts x cells

  if (strip_version) {
    stripped <- sub("\\.[0-9]+$", "", rownames(m))
    n_dup <- sum(duplicated(stripped))
    if (n_dup > 0L) {
      grp <- factor(stripped, levels = unique(stripped))
      agg <- Matrix::sparseMatrix(
        i = as.integer(grp), j = seq_len(nrow(m)), x = 1,
        dims = c(nlevels(grp), nrow(m)))
      m <- agg %*% m
      rownames(m) <- levels(grp)
      rlang::warn(paste0(
        "Stripping transcript version suffixes collapsed ", n_dup,
        " ID(s); summed their counts."))
    } else {
      rownames(m) <- stripped
    }
  }

  if (!is.null(cells)) {
    keep <- intersect(cells, colnames(m))
    if (length(keep) == 0L) {
      rlang::abort(paste0(
        "None of the requested `cells` are present in the matrix barcodes. ",
        "Example requested: ", paste(utils::head(cells, 3L), collapse = ", "),
        "; example matrix: ",
        paste(utils::head(colnames(m), 3L), collapse = ", ")))
    }
    if (length(keep) < length(cells) && verbose) {
      cli::cli_alert_warning(paste0(
        "{length(keep)}/{length(cells)} requested cells found; ",
        "the rest are absent from the matrix and dropped."))
    }
    m <- m[, keep, drop = FALSE]
  }

  out <- methods::as(m, "CsparseMatrix")
  if (verbose) {
    cli::cli_alert_success(paste0(
      "Transcript matrix: {nrow(out)} transcripts x {ncol(out)} cells."))
  }
  out
}
