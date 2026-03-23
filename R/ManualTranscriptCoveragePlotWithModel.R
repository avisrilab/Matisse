#' Plot transcript coverage together with a transcript model.
#'
#' This function combines the baked-in normalized transcript coverage plot with a
#' transcript model panel. The coverage plot always uses the mandatory
#' `plot_value` column produced upstream by `PrepareTranscriptCoverageTracks()`.
#'
#' @param track.df Data frame of transcript coverage tracks. Must contain genomic
#'   coordinates, `plot_value`, and the requested transcript identifier column.
#' @param model List produced by `BuildGeneModel()`. Passed to
#'   `PlotTranscriptModel()`.
#' @param region Optional `GenomicRanges::GRanges` region to plot. If `NULL`,
#'   the region is inferred from `track.df`.
#' @param transcript.col Character(1). Which transcript identifier column to use;
#'   one of `"transcript_name"` or `"transcript_id"`.
#' @param window Numeric(1) or `NA`. Passed to `ManualTranscriptCoveragePlot()`.
#' @param ymax Optional numeric. Passed to `ManualTranscriptCoveragePlot()`.
#' @param split.assays Logical(1). Passed to `ManualTranscriptCoveragePlot()`.
#' @param assay.scale Character(1). Passed to `ManualTranscriptCoveragePlot()`.
#' @param model.height Numeric(1). Relative height of the transcript model panel.
#'
#' @return A combined patchwork plot containing transcript coverage and transcript
#'   model tracks.
#' @export
ManualTranscriptCoveragePlotWithModel <- function(
    track.df,
    model,
    region = NULL,
    transcript.col = c("transcript_name", "transcript_id"),
    window = NA,
    ymax = NULL,
    split.assays = TRUE,
    assay.scale = "separate",
    model.height = 0.4
) {
  
  # ---- Helper checks ----
  .check_data_frame <- function(x, nm) {
    if (!is.data.frame(x)) {
      stop(nm, " must be a data.frame.", call. = FALSE)
    }
  }
  
  .check_named_list <- function(x, nm) {
    if (!is.list(x)) {
      stop(nm, " must be a list.", call. = FALSE)
    }
  }
  
  .check_required_fields <- function(x, fields, nm) {
    miss <- setdiff(fields, names(x))
    if (length(miss) > 0L) {
      stop(
        nm, " is missing required field(s): ",
        paste(miss, collapse = ", "),
        call. = FALSE
      )
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
  
  .check_granges <- function(x, nm) {
    if (!methods::is(x, "GRanges")) {
      stop(nm, " must be a GenomicRanges::GRanges object.", call. = FALSE)
    }
  }
  
  .check_scalar_logical <- function(x, nm) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop(nm, " must be a logical(1).", call. = FALSE)
    }
  }
  
  .check_scalar_numeric_or_na <- function(x, nm) {
    if (!(length(x) == 1L && (is.na(x) || is.numeric(x)))) {
      stop(nm, " must be numeric(1) or NA.", call. = FALSE)
    }
  }
  
  .check_scalar_character <- function(x, nm) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
      stop(nm, " must be a non-empty character(1).", call. = FALSE)
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
  
  .check_named_list(model, "model")
  .check_required_fields(model, c("exons"), "model")
  
  if (!methods::is(model$exons, "GRanges")) {
    stop("model$exons must be a GenomicRanges::GRanges object.", call. = FALSE)
  }
  
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
  
  .check_scalar_numeric_or_na(window, "window")
  
  if (!is.null(ymax)) {
    if (!is.numeric(ymax) || length(ymax) != 1L || is.na(ymax)) {
      stop("ymax must be NULL or numeric(1).", call. = FALSE)
    }
  }
  
  .check_scalar_logical(split.assays, "split.assays")
  .check_scalar_character(assay.scale, "assay.scale")
  
  if (!is.numeric(model.height) || length(model.height) != 1L || is.na(model.height)) {
    stop("model.height must be numeric(1).", call. = FALSE)
  }
  
  # ---- Infer plotting region if needed ----
  if (is.null(region)) {
    chr <- unique(track.df$seqnames)
    if (length(chr) != 1L) {
      stop("track.df must contain exactly one seqname when region is NULL.", call. = FALSE)
    }
    
    region <- GenomicRanges::GRanges(
      seqnames = chr,
      ranges = IRanges::IRanges(
        start = min(track.df$start),
        end = max(track.df$end)
      )
    )
  }
  
  start.pos <- BiocGenerics::start(region)
  end.pos <- BiocGenerics::end(region)
  
  # ---- Coverage panel ----
  p.cov <- ManualTranscriptCoveragePlot(
    track.df = track.df,
    region = region,
    transcript.col = transcript.col,
    window = window,
    ymax = ymax,
    split.assays = split.assays,
    assay.scale = assay.scale,
    annotation = FALSE
  ) &
    ggplot2::scale_x_continuous(
      limits = c(start.pos, end.pos),
      expand = c(0, 0)
    ) &
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  
  # ---- Transcript model panel ----
  p.model <- PlotTranscriptModel(
    model = model,
    region = region,
    transcript.col = transcript.col,
    show.x.axis = TRUE
  )
  
  # ---- Combine panels ----
  p.cov / p.model + patchwork::plot_layout(heights = c(1, model.height))
}