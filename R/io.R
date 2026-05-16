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
