#' Convert transcript coverage tracks into a cut-matrix list for plotting.
#'
#' This function converts plot-ready transcript coverage tracks into a
#' transcript-named list of group x genomic-position matrices. It uses the
#' mandatory normalized plotting column `plot_value` generated upstream by
#' `PrepareTranscriptCoverageTracks()`.
#'
#' @param track.df Data frame of transcript coverage tracks. Must contain genomic
#'   coordinates, grouping information or fields sufficient to create it, the
#'   required `plot_value` column, and the requested transcript identifier column.
#' @param region Optional `GenomicRanges::GRanges` object defining the plotting
#'   region. If `NULL`, the region is inferred from `track.df`.
#' @param transcript.col Character(1). Which transcript identifier column to use;
#'   one of `"transcript_name"` or `"transcript_id"`.
#'
#' @return A list containing a transcript-named cut-matrix list, grouping factor,
#'   and genomic region.
#' @export
TracksToCutmatList <- function(
    track.df,
    region = NULL,
    transcript.col = c("transcript_name", "transcript_id")
) {
  
  # ---- Helper checks ----
  .check_data_frame <- function(x, nm) {
    if (!is.data.frame(x)) {
      stop(nm, " must be a data.frame.", call. = FALSE)
    }
  }
  
  .check_granges <- function(x, nm) {
    if (!methods::is(x, "GRanges")) {
      stop(nm, " must be a GenomicRanges::GRanges object.", call. = FALSE)
    }
  }
  
  .check_required_cols <- function(x, cols, nm) {
    miss <- setdiff(cols, colnames(x))
    if (length(miss) > 0L) {
      stop(
        nm, " is missing required column(s): ",
        paste(miss, collapse = ", "),
        call. = FALSE
      )
    }
  }
  
  # ---- Match args ----
  transcript.col <- match.arg(transcript.col)
  
  # ---- Validate inputs ----
  .check_data_frame(track.df, "track.df")
  
  .check_required_cols(
    track.df,
    c("seqnames", "start", "end", "plot_value", transcript.col),
    "track.df"
  )
  
  if (!is.null(region)) {
    .check_granges(region, "region")
    if (length(region) != 1L) {
      stop("region must be a GenomicRanges::GRanges of length 1.", call. = FALSE)
    }
  }
  
  if (anyNA(track.df$seqnames) || any(!nzchar(as.character(track.df$seqnames)))) {
    stop("track.df$seqnames contains missing or empty values.", call. = FALSE)
  }
  
  if (anyNA(track.df$start)) {
    stop("track.df$start contains missing values.", call. = FALSE)
  }
  
  if (anyNA(track.df$end)) {
    stop("track.df$end contains missing values.", call. = FALSE)
  }
  
  if (!is.numeric(track.df$start)) {
    stop("track.df$start must be numeric.", call. = FALSE)
  }
  
  if (!is.numeric(track.df$end)) {
    stop("track.df$end must be numeric.", call. = FALSE)
  }
  
  if (any(track.df$start < 1)) {
    stop("track.df$start must contain values >= 1.", call. = FALSE)
  }
  
  if (any(track.df$end < track.df$start)) {
    stop("track.df$end must be >= track.df$start for all rows.", call. = FALSE)
  }
  
  if (anyNA(track.df$plot_value)) {
    stop("track.df$plot_value contains missing values.", call. = FALSE)
  }
  
  if (!is.numeric(track.df$plot_value)) {
    stop("track.df$plot_value must be numeric.", call. = FALSE)
  }
  
  if (anyNA(track.df[[transcript.col]]) || any(!nzchar(as.character(track.df[[transcript.col]])))) {
    stop("track.df$", transcript.col, " contains missing or empty values.", call. = FALSE)
  }
  
  if ("group" %in% colnames(track.df)) {
    if (anyNA(track.df$group) || any(!nzchar(as.character(track.df$group)))) {
      stop("track.df$group contains missing or empty values.", call. = FALSE)
    }
  }
  
  # ---- Prepare plotting data ----
  plot.df <- track.df
  
  if (!("group" %in% colnames(plot.df))) {
    plot.df$group <- "All cells"
  }
  
  if (is.null(region)) {
    chr <- unique(plot.df$seqnames)
    if (length(chr) != 1L) {
      stop("track.df must contain exactly one seqname when region is NULL.", call. = FALSE)
    }
    
    region <- GenomicRanges::GRanges(
      seqnames = chr,
      ranges = IRanges::IRanges(
        start = min(plot.df$start),
        end = max(plot.df$end)
      )
    )
  }
  
  start.pos <- BiocGenerics::start(region)
  end.pos <- BiocGenerics::end(region)
  pos <- start.pos:end.pos
  
  groups <- unique(plot.df$group)
  txs <- unique(plot.df[[transcript.col]])
  
  cutmat.list <- vector("list", length(txs))
  names(cutmat.list) <- txs
  
  for (tx in txs) {
    tx.df <- plot.df[plot.df[[transcript.col]] == tx, , drop = FALSE]
    
    mat <- matrix(
      0,
      nrow = length(groups),
      ncol = length(pos),
      dimnames = list(groups, as.character(pos))
    )
    
    for (i in seq_len(nrow(tx.df))) {
      g <- tx.df$group[i]
      s <- max(tx.df$start[i], start.pos)
      e <- min(tx.df$end[i], end.pos)
      v <- tx.df$plot_value[i]
      
      if (e >= s) {
        idx <- (s:e) - start.pos + 1L
        mat[g, idx] <- mat[g, idx] + v
      }
    }
    
    cutmat.list[[tx]] <- mat
  }
  
  obj.groups <- factor(groups, levels = groups)
  names(obj.groups) <- groups
  
  list(
    cutmat.list = cutmat.list,
    obj.groups = obj.groups,
    region = region
  )
}