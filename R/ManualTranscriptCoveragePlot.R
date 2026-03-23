#' Plot transcript coverage tracks using Signac coverage plotting utilities.
#'
#' This function plots transcript coverage tracks using the mandatory normalized
#' plotting column `plot_value` generated upstream by
#' `PrepareTranscriptCoverageTracks()`. Group scaling and global scaling are
#' therefore assumed to already be baked into `track.df`. Smoothing is applied
#' here through the `window` argument passed to `Signac:::CoverageTrack()`.
#'
#' @param track.df Data frame of transcript coverage tracks. Must contain genomic
#'   coordinates, the required `plot_value` column, and the requested transcript
#'   identifier column.
#' @param region Optional `GenomicRanges::GRanges` region to plot. If `NULL`,
#'   the region is inferred from `track.df`.
#' @param transcript.col Character(1). Which transcript identifier column to use;
#'   one of `"transcript_name"` or `"transcript_id"`.
#' @param window Numeric(1) or `NA`. Passed to `Signac:::CoverageTrack()` for
#'   smoothing.
#' @param ymax Optional numeric. Passed to `Signac:::CoverageTrack()`.
#' @param split.assays Logical(1). Passed to `Signac:::CoverageTrack()`.
#' @param assay.scale Character(1). Passed to `Signac:::CoverageTrack()`.
#' @param annotation Logical(1). If `TRUE`, add a gene annotation track.
#' @param object Optional Seurat object. Required if `annotation = TRUE`.
#' @param assay Optional Character(1). Assay name to use for annotation plotting.
#'   Required if `annotation = TRUE`.
#' @param heights Optional numeric vector. Passed to `Signac:::CombineTracks()`.
#'
#' @return A ggplot/patchwork-style combined coverage plot.
#' @export
ManualTranscriptCoveragePlot <- function(
    track.df,
    region = NULL,
    transcript.col = c("transcript_name", "transcript_id"),
    window = NA,
    ymax = NULL,
    split.assays = TRUE,
    assay.scale = "separate",
    annotation = FALSE,
    object = NULL,
    assay = NULL,
    heights = NULL
) {
  
  # ---- Helper checks ----
  .check_data_frame <- function(x, nm) {
    if (!is.data.frame(x)) {
      stop(nm, " must be a data.frame.", call. = FALSE)
    }
  }
  
  .check_scalar_logical <- function(x, nm) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop(nm, " must be a logical(1).", call. = FALSE)
    }
  }
  
  .check_scalar_character <- function(x, nm) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
      stop(nm, " must be a non-empty character(1).", call. = FALSE)
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
  
  if (!(length(window) == 1L && (is.na(window) || is.numeric(window)))) {
    stop("window must be numeric(1) or NA.", call. = FALSE)
  }
  
  if (!is.null(ymax)) {
    if (!is.numeric(ymax) || length(ymax) != 1L || is.na(ymax)) {
      stop("ymax must be NULL or numeric(1).", call. = FALSE)
    }
  }
  
  .check_scalar_logical(split.assays, "split.assays")
  .check_scalar_character(assay.scale, "assay.scale")
  .check_scalar_logical(annotation, "annotation")
  
  if (!is.null(heights)) {
    if (!is.numeric(heights) || anyNA(heights)) {
      stop("heights must be NULL or a numeric vector with no missing values.", call. = FALSE)
    }
  }
  
  if (isTRUE(annotation)) {
    if (is.null(object)) {
      stop("If annotation = TRUE, object must be provided.", call. = FALSE)
    }
    if (!isS4(object)) {
      stop("object must be a Seurat object.", call. = FALSE)
    }
    if (is.null(assay)) {
      stop("If annotation = TRUE, assay must be provided.", call. = FALSE)
    }
    .check_scalar_character(assay, "assay")
    assay.names <- names(object@assays)
    if (!(assay %in% assay.names)) {
      stop("assay not found in object: ", assay, call. = FALSE)
    }
  } else {
    if (!is.null(assay)) {
      .check_scalar_character(assay, "assay")
    }
    if (!is.null(object) && !isS4(object)) {
      stop("object must be NULL or a Seurat object.", call. = FALSE)
    }
  }
  
  # ---- Build cut matrices from baked-in normalized plotting values ----
  cm.obj <- TracksToCutmatList(
    track.df = track.df,
    region = region,
    transcript.col = transcript.col
  )
  
  cutmat.list <- cm.obj$cutmat.list
  obj.groups <- cm.obj$obj.groups
  region <- cm.obj$region
  
  # ---- Scaling already applied upstream; pass neutral factors to CoverageTrack ----
  gsf.list <- lapply(cutmat.list, function(x) {
    out <- rep(1, length(levels(obj.groups)))
    names(out) <- levels(obj.groups)
    out
  })
  
  sf.list <- lapply(cutmat.list, function(x) 1)
  
  p <- Signac:::CoverageTrack(
    cutmat = cutmat.list,
    region = region,
    group.scale.factors = gsf.list,
    scale.factor = sf.list,
    window = window,
    ymax = ymax,
    split.assays = split.assays,
    assay.scale = assay.scale,
    obj.groups = obj.groups,
    downsample.rate = 1,
    max.downsample = Inf
  )
  
  if (isTRUE(annotation)) {
    gene.plot <- Signac::AnnotationPlot(
      object = object[[assay]],
      region = region,
      mode = "gene"
    )
  } else {
    gene.plot <- NULL
  }
  
  heights <- Signac:::SetIfNull(
    x = heights,
    y = c(10, 3)
  )
  
  p <- Signac:::CombineTracks(
    plotlist = list(p, gene.plot),
    heights = heights
  ) & ggplot2::theme(
    legend.key.size = grid::unit(1 / 2, "lines"),
    legend.text = ggplot2::element_text(size = 7),
    legend.title = ggplot2::element_text(size = 8)
  )
  
  return(p)
}