#' Plot a transcript model from exon structures.
#'
#' @param model List produced by `BuildGeneModel()`. Must contain `exons` with
#'   transcript identifier metadata columns `tx_id` and `tx_name`.
#' @param region Optional `GenomicRanges::GRanges` region to plot. If `NULL`,
#'   the plotting range is inferred from the exon coordinates.
#' @param transcript.col Character(1). Which transcript label to display; one of
#'   `"transcript_name"` or `"transcript_id"`.
#' @param show.x.axis Logical(1). Whether to display the x-axis.
#'
#' @return A ggplot object showing exon blocks and intron connectors for each
#'   transcript.
#' @export
PlotTranscriptModel <- function(
    model,
    region = NULL,
    transcript.col = c("transcript_name", "transcript_id"),
    show.x.axis = TRUE
) {
  
  # ---- Helper checks ----
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
  
  # ---- Match args ----
  transcript.col <- match.arg(transcript.col)
  
  # ---- Validate inputs ----
  .check_named_list(model, "model")
  .check_required_fields(model, c("exons"), "model")
  
  .check_granges(model$exons, "model$exons")
  .check_scalar_logical(show.x.axis, "show.x.axis")
  
  if (!is.null(region)) {
    .check_granges(region, "region")
    if (length(region) != 1L) {
      stop("region must be a GenomicRanges::GRanges of length 1.", call. = FALSE)
    }
  }
  
  if (length(model$exons) < 1L) {
    stop("model$exons must contain at least one range.", call. = FALSE)
  }
  
  exon.mcols <- S4Vectors::mcols(model$exons)
  required.exon.cols <- c("tx_id", "tx_name")
  miss.exon.cols <- setdiff(required.exon.cols, colnames(exon.mcols))
  if (length(miss.exon.cols) > 0L) {
    stop(
      "model$exons is missing required metadata column(s): ",
      paste(miss.exon.cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (anyNA(exon.mcols$tx_id) || any(!nzchar(as.character(exon.mcols$tx_id)))) {
    stop("model$exons$tx_id contains missing or empty values.", call. = FALSE)
  }
  
  if (anyNA(exon.mcols$tx_name) || any(!nzchar(as.character(exon.mcols$tx_name)))) {
    stop("model$exons$tx_name contains missing or empty values.", call. = FALSE)
  }
  
  # ---- Build exon data frame ----
  exons <- model$exons
  
  exon.df <- data.frame(
    transcript_id = S4Vectors::mcols(exons)$tx_id,
    transcript_name = S4Vectors::mcols(exons)$tx_name,
    seqnames = as.character(GenomicRanges::seqnames(exons)),
    start = BiocGenerics::start(exons),
    end = BiocGenerics::end(exons),
    strand = as.character(BiocGenerics::strand(exons)),
    stringsAsFactors = FALSE
  )
  
  if (is.null(region)) {
    start.pos <- min(exon.df$start)
    end.pos <- max(exon.df$end)
  } else {
    start.pos <- BiocGenerics::start(region)
    end.pos <- BiocGenerics::end(region)
  }
  
  exon.df$label <- exon.df[[transcript.col]]
  tx.order <- rev(unique(exon.df$label))
  exon.df$label <- factor(exon.df$label, levels = tx.order)
  
  # ---- Build intron data frame ----
  exon.list <- split(exons, S4Vectors::mcols(exons)$tx_id)
  
  intron.list <- lapply(exon.list, function(gr) {
    gr <- sort(gr)
    if (length(gr) < 2L) {
      return(NULL)
    }
    
    data.frame(
      transcript_id = unique(S4Vectors::mcols(gr)$tx_id),
      transcript_name = unique(S4Vectors::mcols(gr)$tx_name),
      start = BiocGenerics::end(gr)[-length(gr)],
      end = BiocGenerics::start(gr)[-1],
      stringsAsFactors = FALSE
    )
  })
  
  keep.introns <- !vapply(intron.list, is.null, logical(1))
  if (any(keep.introns)) {
    intron.df <- do.call(rbind, intron.list[keep.introns])
    rownames(intron.df) <- NULL
  } else {
    intron.df <- data.frame(
      transcript_id = character(0),
      transcript_name = character(0),
      start = numeric(0),
      end = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  intron.df$label <- intron.df[[transcript.col]]
  intron.df$label <- factor(intron.df$label, levels = tx.order)
  
  y.map <- seq_along(tx.order)
  names(y.map) <- tx.order
  
  exon.df$y <- unname(y.map[as.character(exon.df$label)])
  intron.df$y <- unname(y.map[as.character(intron.df$label)])
  
  # ---- Plot ----
  p <- ggplot2::ggplot()
  
  if (nrow(intron.df) > 0L) {
    p <- p +
      ggplot2::geom_segment(
        data = intron.df,
        ggplot2::aes(
          x = start,
          xend = end,
          y = y,
          yend = y
        ),
        linewidth = 0.3,
        lineend = "round"
      )
  }
  
  p <- p +
    ggplot2::geom_rect(
      data = exon.df,
      ggplot2::aes(
        xmin = start,
        xmax = end,
        ymin = y - 0.25,
        ymax = y + 0.25
      ),
      fill = "black"
    ) +
    ggplot2::scale_y_continuous(
      breaks = y.map,
      labels = names(y.map),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::scale_x_continuous(
      limits = c(start.pos, end.pos),
      expand = c(0, 0)
    ) +
    ggplot2::labs(x = "Genomic position", y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 8)
    )
  
  if (!show.x.axis) {
    p <- p +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }
  
  p
}