#' Prepare transcript coverage tracks for plotting.
#'
#' This function aggregates transcript-resolved signal and read support to the
#' group x transcript x bin level, computes mandatory Signac-style normalization,
#' and returns a plot-ready track table with a single `plot_value` column used
#' downstream for visualization.
#'
#' @param signal.df Data frame of transcript-resolved signal. Must contain at least
#'   `value`, `transcript_id`, `transcript_name`, and `bin_id`. If `group` is not
#'   present, it will be added as `"All cells"`.
#' @param read.df Data frame of transcript-resolved read support. Must contain at
#'   least `value`, `transcript_id`, `transcript_name`, and `bin_id`. If `group`
#'   is not present, it will be added as `"All cells"`.
#' @param model List produced by `BuildGeneModel()`. Must contain `bins` with
#'   required bin-level metadata.
#'
#' @return A data frame containing aggregated signal, read support, mandatory
#'   normalization factors, and genomic bin coordinates for plotting. The column
#'   `plot_value` is the only value intended for downstream plotting.
#' @export
PrepareTranscriptCoverageTracks <- function(
    signal.df,
    read.df,
    model
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
  
  .check_support_df <- function(x, nm) {
    .check_data_frame(x, nm)
    .check_required_fields(
      x,
      c("value", "transcript_id", "transcript_name", "bin_id"),
      nm
    )
    
    if (anyNA(x$value)) {
      stop(nm, "$value contains missing values.", call. = FALSE)
    }
    if (!is.numeric(x$value)) {
      stop(nm, "$value must be numeric.", call. = FALSE)
    }
    if (anyNA(x$transcript_id) || any(!nzchar(as.character(x$transcript_id)))) {
      stop(nm, "$transcript_id contains missing or empty values.", call. = FALSE)
    }
    if (anyNA(x$transcript_name) || any(!nzchar(as.character(x$transcript_name)))) {
      stop(nm, "$transcript_name contains missing or empty values.", call. = FALSE)
    }
    if (anyNA(x$bin_id) || any(!nzchar(as.character(x$bin_id)))) {
      stop(nm, "$bin_id contains missing or empty values.", call. = FALSE)
    }
    if ("group" %in% colnames(x)) {
      if (anyNA(x$group) || any(!nzchar(as.character(x$group)))) {
        stop(nm, "$group contains missing or empty values.", call. = FALSE)
      }
    }
  }
  
  # ---- Validate inputs ----
  .check_support_df(signal.df, "signal.df")
  .check_support_df(read.df, "read.df")
  
  .check_named_list(model, "model")
  .check_required_fields(model, c("bins"), "model")
  
  if (!methods::is(model$bins, "GRanges")) {
    stop("model$bins must be a GenomicRanges::GRanges object.", call. = FALSE)
  }
  
  if (length(model$bins) < 1L) {
    stop("model$bins must contain at least one range.", call. = FALSE)
  }
  
  bin.mcols <- S4Vectors::mcols(model$bins)
  required.bin.cols <- c("bin_id", "n_tx", "is_shared")
  miss.bin.cols <- setdiff(required.bin.cols, colnames(bin.mcols))
  if (length(miss.bin.cols) > 0L) {
    stop(
      "model$bins is missing required metadata column(s): ",
      paste(miss.bin.cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (anyNA(bin.mcols$bin_id) || any(!nzchar(as.character(bin.mcols$bin_id)))) {
    stop("model$bins$bin_id contains missing or empty values.", call. = FALSE)
  }
  
  if (anyNA(bin.mcols$n_tx)) {
    stop("model$bins$n_tx contains missing values.", call. = FALSE)
  }
  
  if (anyNA(bin.mcols$is_shared)) {
    stop("model$bins$is_shared contains missing values.", call. = FALSE)
  }
  
  # ---- add default group if missing ----
  if (!("group" %in% colnames(signal.df))) {
    signal.df$group <- "All cells"
  }
  if (!("group" %in% colnames(read.df))) {
    read.df$group <- "All cells"
  }
  
  # ---- keep only rows with positive support ----
  signal.df <- signal.df[signal.df$value > 0, , drop = FALSE]
  read.df <- read.df[read.df$value > 0, , drop = FALSE]
  
  # ---- aggregate signal per group x transcript x bin ----
  signal.agg <- stats::aggregate(
    value ~ group + transcript_id + transcript_name + bin_id,
    data = signal.df,
    FUN = sum
  )
  colnames(signal.agg)[colnames(signal.agg) == "value"] <- "signal_value"
  
  # ---- aggregate read support per group x transcript x bin ----
  read.agg <- stats::aggregate(
    value ~ group + transcript_id + transcript_name + bin_id,
    data = read.df,
    FUN = sum
  )
  colnames(read.agg)[colnames(read.agg) == "value"] <- "read_support"
  
  # ---- merge aggregated signal and read support ----
  agg.df <- merge(
    signal.agg,
    read.agg,
    by = c("group", "transcript_id", "transcript_name", "bin_id"),
    all = TRUE,
    sort = FALSE
  )
  
  agg.df$signal_value[is.na(agg.df$signal_value)] <- 0
  agg.df$read_support[is.na(agg.df$read_support)] <- 0
  
  # ---- compute group-level scale factors from total read support ----
  group.depth <- stats::aggregate(
    read_support ~ group,
    data = agg.df,
    FUN = sum
  )
  colnames(group.depth)[colnames(group.depth) == "read_support"] <- "group_scale_factor"
  
  # ---- choose global reference scale ----
  ref.scale <- stats::median(group.depth$group_scale_factor[group.depth$group_scale_factor > 0])
  if (length(ref.scale) == 0L || is.na(ref.scale)) {
    ref.scale <- 1
  }
  
  group.depth$global_scale_factor <- ref.scale
  
  # ---- merge scale factors back ----
  agg.df <- merge(
    agg.df,
    group.depth,
    by = "group",
    all.x = TRUE,
    sort = FALSE
  )
  
  # ---- mandatory normalized plotting value ----
  agg.df$plot_value <- ifelse(
    agg.df$group_scale_factor > 0,
    agg.df$signal_value * (agg.df$global_scale_factor / agg.df$group_scale_factor),
    0
  )
  
  # ---- retain raw signal for bookkeeping ----
  agg.df$raw_value <- agg.df$signal_value
  
  # ---- add bin genomic coordinates ----
  bin.df <- data.frame(
    bin_id = S4Vectors::mcols(model$bins)$bin_id,
    seqnames = as.character(GenomicRanges::seqnames(model$bins)),
    start = BiocGenerics::start(model$bins),
    end = BiocGenerics::end(model$bins),
    strand = as.character(BiocGenerics::strand(model$bins)),
    n_tx = S4Vectors::mcols(model$bins)$n_tx,
    is_shared = S4Vectors::mcols(model$bins)$is_shared,
    stringsAsFactors = FALSE
  )
  
  track.df <- merge(
    agg.df,
    bin.df,
    by = "bin_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  # ---- order rows for plotting ----
  track.df <- track.df[order(
    track.df$group,
    track.df$transcript_name,
    track.df$start,
    track.df$end
  ), , drop = FALSE]
  
  rownames(track.df) <- NULL
  
  track.df
}