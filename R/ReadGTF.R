#' Read a GTF and extract transcript-level annotations aligned to a transcript matrix
#'
#' @param tmat A sparse matrix with transcripts/isoforms as rows and cells as columns.
#' @param gtf.file Character(1). Path to a GTF file (.gtf or .gtf.gz).
#' @param transcript.field Character(1). Attribute key in the GTF to match transcripts in `tmat` (default = "transcript_id").
#' @param gene.field Character(1). Attribute key for the gene label to return (default = "gene_name").
#' @param transcript.type.field Character(1). Attribute key for transcript type/biotype to return (default = "transcript_type").
#' @param include.gene_id Logical(1). If TRUE, also extracts `gene_id` into a `gene.id` column (default = FALSE).
#' @param verbose Logical(1). If TRUE, prints basic progress messages (default = FALSE).
#'
#' @return A data.frame with one row per transcript in `tmat`, containing:
#'   `transcript_id`, `gene`, `transcript.type`, `seqname`, `start`, `end`, `strand`, and optional `gene.id`. 
#' @export
ReadGTF <- function(tmat,
                    gtf.file,
                    transcript.field = "transcript_id",
                    gene.field = "gene_name",
                    transcript.type.field = "transcript_type",
                    include.gene_id = FALSE,
                    verbose = FALSE) {
  
  # ---- Helper functions ----
  .check_scalar_string <- function(x, nm) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
      stop(nm, " must be a non-empty character(1).", call. = FALSE)
    }
  }
  
  .check_flag <- function(x, nm) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop(nm, " must be TRUE or FALSE.", call. = FALSE)
    }
  }
  
  .extract_attr <- function(x, key) {
    pat <- paste0('(?:^|;[[:space:]]*)', key, '[[:space:]]+"([^"]+)"')
    m <- regexec(pat, x, perl = TRUE)
    regm <- regmatches(x, m)
    vapply(regm, function(z) {
      if (length(z) >= 2L) z[2] else NA_character_
    }, character(1))
  }
  
  .open_text <- function(path) {
    if (grepl("\\.gz$", path, ignore.case = TRUE)) {
      gzfile(path, open = "rt")
    } else {
      file(path, open = "rt")
    }
  }
  
  # ---- Verify parameter inputs ----
  if (!(inherits(tmat, "Matrix") || is.matrix(tmat))) {
    stop("tmat must be a matrix-like object with transcripts as rows and cells as columns.",
         call. = FALSE)
  }
  
  tx.in.mat <- rownames(tmat)
  
  if (is.null(tx.in.mat) || anyNA(tx.in.mat) || any(!nzchar(tx.in.mat))) {
    stop("tmat must have non-empty rownames containing transcript identifiers.",
         call. = FALSE)
  }
  
  tx.in.mat <- trimws(as.character(tx.in.mat))
  
  if (anyNA(tx.in.mat) || any(!nzchar(tx.in.mat))) {
    stop("tmat rownames contain NA or empty transcript identifiers after trimming.",
         call. = FALSE)
  }
  
  if (anyDuplicated(tx.in.mat)) {
    warning("Duplicate transcript identifiers detected in tmat rownames.",
            call. = FALSE)
  }
  
  .check_scalar_string(gtf.file, "gtf.file")
  .check_scalar_string(transcript.field, "transcript.field")
  .check_scalar_string(gene.field, "gene.field")
  .check_scalar_string(transcript.type.field, "transcript.type.field")
  .check_flag(include.gene_id, "include.gene_id")
  .check_flag(verbose, "verbose")
  
  if (!file.exists(gtf.file)) {
    stop("File not found: ", gtf.file, call. = FALSE)
  }
  
  if (!grepl("\\.gtf(?:\\.gz)?$", gtf.file, ignore.case = TRUE)) {
    warning("gtf.file does not end in .gtf or .gtf.gz; attempting to parse as GTF anyway.",
            call. = FALSE)
  }
  
  # ---- Read GTF ----
  if (isTRUE(verbose)) {
    message("Reading GTF: ", gtf.file)
  }
  
  con <- .open_text(gtf.file)
  on.exit(close(con), add = TRUE)
  
  gtf <- tryCatch(
    utils::read.delim(
      con,
      header = FALSE,
      sep = "\t",
      quote = "",
      comment.char = "#",
      stringsAsFactors = FALSE
    ),
    error = function(e) {
      stop("Failed to read GTF: ", conditionMessage(e), call. = FALSE)
    }
  )
  
  if (nrow(gtf) < 1L) {
    stop("GTF contains no data rows.", call. = FALSE)
  }
  
  if (ncol(gtf) != 9L) {
    stop("Invalid GTF: expected exactly 9 tab-separated columns, found ",
         ncol(gtf), ".", call. = FALSE)
  }
  
  colnames(gtf) <- c(
    "seqname", "source", "feature", "start", "end",
    "score", "strand", "frame", "attributes"
  )
  
  if (!all(nzchar(as.character(gtf$seqname)))) {
    stop("Invalid GTF: seqname column contains empty values.", call. = FALSE)
  }
  
  gtf$start <- suppressWarnings(as.integer(gtf$start))
  gtf$end   <- suppressWarnings(as.integer(gtf$end))
  
  if (anyNA(gtf$start) || anyNA(gtf$end)) {
    stop("Invalid GTF: start/end coordinates could not be parsed as integers.",
         call. = FALSE)
  }
  
  if (any(gtf$start < 1L) || any(gtf$end < 1L)) {
    stop("Invalid GTF: genomic coordinates must be positive integers.",
         call. = FALSE)
  }
  
  if (any(gtf$end < gtf$start)) {
    stop("Invalid GTF: found rows where end < start.", call. = FALSE)
  }
  
  valid.strand <- c("+", "-", ".")
  if (any(!gtf$strand %in% valid.strand)) {
    stop("Invalid GTF: strand column must contain only '+', '-', or '.'.",
         call. = FALSE)
  }
  
  if (anyNA(gtf$attributes) || any(!nzchar(gtf$attributes))) {
    stop("Invalid GTF: attributes column contains missing or empty entries.",
         call. = FALSE)
  }
  
  idx <- which(gtf$feature == "transcript")
  if (length(idx) == 0L) {
    stop("Invalid GTF: no 'transcript' feature rows found.", call. = FALSE)
  }
  
  sub <- gtf[idx, c("seqname", "start", "end", "strand", "attributes"), drop = FALSE]
  
  tx <- trimws(.extract_attr(sub$attributes, transcript.field))
  gn <- trimws(.extract_attr(sub$attributes, gene.field))
  tt <- trimws(.extract_attr(sub$attributes, transcript.type.field))
  
  if (all(is.na(tx) | !nzchar(tx))) {
    stop("Invalid GTF: no transcript rows contain the requested transcript field '",
         transcript.field, "'.", call. = FALSE)
  }
  
  if (all(is.na(tt) | !nzchar(tt))) {
    alt.key <- "transcript_biotype"
    if (!identical(transcript.type.field, alt.key)) {
      tt.alt <- trimws(.extract_attr(sub$attributes, alt.key))
      if (!all(is.na(tt.alt) | !nzchar(tt.alt))) {
        warning("'", transcript.type.field, "' not found in GTF attributes; using '",
                alt.key, "' instead.", call. = FALSE)
        tt <- tt.alt
      } else {
        warning("Transcript type field not found in GTF attributes (tried '",
                transcript.type.field, "' and '", alt.key,
                "'). Returning transcript.type as NA.", call. = FALSE)
        tt <- rep(NA_character_, length(tx))
      }
    } else {
      warning("Transcript type field '", transcript.type.field,
              "' not found in GTF attributes. Returning transcript.type as NA.",
              call. = FALSE)
      tt <- rep(NA_character_, length(tx))
    }
  }
  
  # ---- Create GTF DF ----
  gtf.df <- data.frame(
    transcript_id = tx,
    gene = gn,
    transcript.type = tt,
    seqname = as.character(sub$seqname),
    start = sub$start,
    end = sub$end,
    strand = as.character(sub$strand),
    stringsAsFactors = FALSE
  )
  
  if (isTRUE(include.gene_id)) {
    gtf.df$gene.id <- trimws(.extract_attr(sub$attributes, "gene_id"))
  }
  
  keep <- !is.na(gtf.df$transcript_id) & nzchar(gtf.df$transcript_id)
  n.drop <- sum(!keep)
  
  if (n.drop > 0L) {
    warning(n.drop, " transcript rows were dropped because transcript_id could not be parsed.",
            call. = FALSE)
  }
  
  gtf.df <- gtf.df[keep, , drop = FALSE]
  
  if (nrow(gtf.df) < 1L) {
    stop("No valid transcript annotations remained after parsing the GTF.",
         call. = FALSE)
  }
  
  if (anyDuplicated(gtf.df$transcript_id)) {
    warning("Duplicate transcript IDs detected in GTF; keeping first occurrence.",
            call. = FALSE)
    gtf.df <- gtf.df[!duplicated(gtf.df$transcript_id), , drop = FALSE]
  }
  
  match.idx <- match(tx.in.mat, gtf.df$transcript_id)
  out <- gtf.df[match.idx, , drop = FALSE]
  
  out$transcript_id <- tx.in.mat
  
  n.missing <- sum(is.na(match.idx))
  if (n.missing > 0L) {
    warning(n.missing,
            " transcripts in tmat were not found in the GTF and were retained with NA annotations.",
            call. = FALSE)
  }
  
  rownames(out) <- NULL
  
  if (isTRUE(verbose)) {
    message("Returning ", nrow(out), " transcript annotations aligned to tmat.")
  }
  
  out
}