#' @include MatisseObject-class.R
#' @include generics.R
NULL

# Shared ggplot2 theme for Matisse plots
.matisse_theme <- function() {
  ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40"),
      legend.position = "right"
    )
}

# ---------------------------------------------------------------------------
# PlotPSIUMAP
# ---------------------------------------------------------------------------

#' UMAP plot coloured by PSI of a splice event
#'
#' Overlays the PSI value of a single splice event on the UMAP embedding
#' stored in the embedded Seurat object.
#'
#' @param object A \code{MatisseObject} with a \code{"psi"} assay.
#'   The Seurat object must contain a \code{umap} reduction.
#' @param event_id Character. Column name in the PSI matrix.
#' @param reduction Character. Name of the dimensionality reduction to use.
#'   Default: \code{"umap"}.
#' @param dims Integer vector of length 2 selecting which dimensions to plot.
#'   Default: \code{c(1, 2)}.
#' @param pt_size Numeric. Point size. Default: \code{0.5}.
#' @param na_colour Character. Colour for \code{NA} PSI cells.
#'   Default: \code{"grey80"}.
#' @param title Character. Plot title. Defaults to the event ID.
#'
#' @return A \code{ggplot} object.
#'
#' @rdname PlotPSIUMAP
#' @export
setMethod("PlotPSIUMAP", "MatisseObject",
          function(object, event_id,
                   reduction  = "umap",
                   dims       = c(1L, 2L),
                   pt_size    = 0.5,
                   na_colour  = "grey80",
                   title      = NULL) {
  .require_psi(object)
  .require_event(object, event_id)

  emb       <- SeuratObject::Embeddings(object@seurat, reduction = reduction)
  dim_names <- colnames(emb)[dims]
  psi_cx    <- GetPSI(object)  # cells x events
  psi_vals  <- as.numeric(psi_cx[rownames(emb), event_id])

  df <- data.frame(
    x   = emb[, dims[1]],
    y   = emb[, dims[2]],
    psi = psi_vals
  )

  plot_title <- title %||% event_id

  ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y,
                                    colour = .data$psi)) +
    ggplot2::geom_point(size = pt_size, na.rm = TRUE) +
    ggplot2::scale_colour_gradientn(
      colours  = c("#2166ac", "#f7f7f7", "#d6604d"),
      na.value = na_colour,
      limits   = c(0, 1),
      name     = "PSI"
    ) +
    ggplot2::labs(
      title = plot_title,
      x     = dim_names[1],
      y     = dim_names[2]
    ) +
    .matisse_theme()
})

# ---------------------------------------------------------------------------
# PlotPSIViolin
# ---------------------------------------------------------------------------

#' Violin plot of PSI values split by cell group
#'
#' @param object A \code{MatisseObject} with a \code{"psi"} assay.
#' @param event_id Character. Column name in the PSI matrix.
#' @param group_by Character. Column in \code{Seurat::meta.data} to split
#'   cells by. Default: \code{"seurat_clusters"}.
#' @param colours Named character vector mapping group levels to colours.
#'   Default: \code{NULL} (uses ggplot2 defaults).
#' @param add_points Logical. Overlay individual cell PSI values as jittered
#'   points. Default: \code{FALSE}.
#' @param title Character. Plot title. Defaults to the event ID.
#'
#' @return A \code{ggplot} object.
#'
#' @rdname PlotPSIViolin
#' @export
setMethod("PlotPSIViolin", "MatisseObject",
          function(object, event_id,
                   group_by   = "seurat_clusters",
                   colours    = NULL,
                   add_points = FALSE,
                   title      = NULL) {
  .require_psi(object)
  .require_event(object, event_id)

  cells    <- .get_cells(object)
  psi_cx   <- GetPSI(object)
  psi_vals <- as.numeric(psi_cx[cells, event_id])
  group_vals <- .get_seurat_meta_col(object, group_by)

  df <- data.frame(psi = psi_vals, group = group_vals)
  df <- df[!is.na(df$psi), ]

  plot_title <- title %||% event_id

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data$group, y = .data$psi, fill = .data$group)) +
    ggplot2::geom_violin(trim = FALSE, scale = "width") +
    ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
    ggplot2::scale_y_continuous(limits = c(0, 1), name = "PSI") +
    ggplot2::labs(title = plot_title, x = group_by) +
    ggplot2::theme(legend.position = "none") +
    .matisse_theme()

  if (!is.null(colours)) p <- p + ggplot2::scale_fill_manual(values = colours)

  if (add_points) {
    p <- p + ggplot2::geom_jitter(
      width = 0.15, size = 0.3, alpha = 0.4,
      ggplot2::aes(colour = .data$group), show.legend = FALSE)
  }

  p
})

# ---------------------------------------------------------------------------
# PlotPSIHeatmap
# ---------------------------------------------------------------------------

#' Heatmap of PSI values (cells x events)
#'
#' @param object A \code{MatisseObject} with a \code{"psi"} assay.
#' @param events Character vector of event IDs to include.
#'   Default: all events.
#' @param cells Character vector of cell barcodes to include.
#'   Default: all cells.
#' @param group_by Character. Column in \code{Seurat::meta.data} used to
#'   annotate and order cells. Default: \code{NULL}.
#' @param max_cells Integer. Downsample to this many cells before plotting.
#'   Default: \code{500}.
#' @param na_colour Character. Colour for \code{NA} entries.
#'   Default: \code{"grey90"}.
#'
#' @return A \code{ggplot} object.
#'
#' @rdname PlotPSIHeatmap
#' @export
setMethod("PlotPSIHeatmap", "MatisseObject",
          function(object,
                   events     = NULL,
                   cells      = NULL,
                   group_by   = NULL,
                   max_cells  = 500L,
                   na_colour  = "grey90") {
  .require_psi(object)

  psi_cx     <- GetPSI(object)  # cells x events
  all_cells  <- .get_cells(object)
  all_events <- colnames(psi_cx)

  cells  <- cells  %||% all_cells
  events <- events %||% all_events

  bad_cells  <- setdiff(cells,  all_cells)
  bad_events <- setdiff(events, all_events)
  if (length(bad_cells)  > 0) rlang::warn("Some cell barcodes not found; skipping.")
  if (length(bad_events) > 0) rlang::warn("Some event IDs not found; skipping.")
  cells  <- intersect(cells,  all_cells)
  events <- intersect(events, all_events)

  if (length(cells) > max_cells) cells <- sample(cells, max_cells)

  psi_sub <- as.matrix(psi_cx[cells, events, drop = FALSE])

  finite_mask <- is.finite(psi_sub)
  if (sum(finite_mask) > 0) {
    psi_for_clust           <- psi_sub
    psi_for_clust[!finite_mask] <- 0.5
    event_order <- tryCatch({
      hc <- stats::hclust(stats::dist(t(psi_for_clust)))
      hc$order
    }, error = function(e) seq_len(ncol(psi_sub)))
  } else {
    event_order <- seq_len(ncol(psi_sub))
  }

  if (!is.null(group_by)) {
    grp   <- .get_seurat_meta_col(object, group_by)[match(cells, all_cells)]
    cells <- cells[order(grp)]
  }

  df <- data.frame(
    cell  = rep(cells,                 each  = length(events)),
    event = rep(events[event_order],   times = length(cells)),
    psi   = as.vector(psi_sub[cells, events[event_order]])
  )
  df$cell  <- factor(df$cell,  levels = cells)
  df$event <- factor(df$event, levels = events[event_order])

  ggplot2::ggplot(df, ggplot2::aes(
    x    = .data$event,
    y    = .data$cell,
    fill = .data$psi)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colours  = c("#2166ac", "#f7f7f7", "#d6604d"),
      na.value = na_colour,
      limits   = c(0, 1),
      name     = "PSI"
    ) +
    ggplot2::labs(title = "PSI Heatmap", x = "Splice Event", y = "Cell") +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    ) +
    .matisse_theme()
})

# ---------------------------------------------------------------------------
# PlotJunctionCoverage
# ---------------------------------------------------------------------------

#' Junction coverage bar plot for a gene
#'
#' @param object A \code{MatisseObject} with non-\code{NULL}
#'   \code{junction_counts} and \code{junction_data} slots.
#' @param gene Character. Gene name to plot.
#' @param cells Character vector of cell barcodes to aggregate over.
#'   Default: all cells.
#' @param log_scale Logical. Use log10 y-axis. Default: \code{FALSE}.
#'
#' @return A \code{ggplot} object.
#'
#' @rdname PlotJunctionCoverage
#' @export
setMethod("PlotJunctionCoverage", "MatisseObject",
          function(object, gene, cells = NULL, log_scale = FALSE) {
  if (is.null(object@junction_counts)) {
    rlang::abort("junction_counts slot is NULL.")
  }
  if (nrow(object@junction_data) == 0) {
    rlang::abort("junction_data is empty.")
  }

  jd        <- object@junction_data
  gene_jxns <- jd$junction_id[jd$gene_id == gene]
  if (length(gene_jxns) == 0) {
    rlang::abort(paste0("No junctions found for gene: ", gene))
  }

  present <- intersect(gene_jxns, colnames(object@junction_counts))
  if (length(present) == 0) {
    rlang::abort(paste0("No junctions for '", gene,
                        "' are present in junction_counts."))
  }

  sub_cells <- cells %||% .get_cells(object)
  sub_mat   <- object@junction_counts[sub_cells, present, drop = FALSE]
  totals    <- Matrix::colSums(sub_mat)

  df         <- data.frame(junction = names(totals), count = as.numeric(totals))
  coord_info <- jd[match(df$junction, jd$junction_id), ]
  df$start   <- coord_info$start
  df         <- df[order(df$start), ]
  df$junction <- factor(df$junction, levels = df$junction)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$junction, y = .data$count)) +
    ggplot2::geom_col(fill = "#4393c3") +
    ggplot2::labs(
      title = paste0("Junction coverage: ", gene),
      x     = "Junction (ordered by genomic position)",
      y     = "Total reads"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7)
    ) +
    .matisse_theme()

  if (log_scale) p <- p + ggplot2::scale_y_log10()
  p
})

# ---------------------------------------------------------------------------
# PlotQCMetrics
# ---------------------------------------------------------------------------

#' Violin plots of isoform QC metrics
#'
#' @param object A \code{MatisseObject} with populated \code{isoform_metadata}.
#' @param features Character vector of QC metric names. Must be columns in
#'   \code{isoform_metadata}. Default: all numeric columns.
#' @param group_by Character. Column in \code{Seurat::meta.data} to split
#'   cells by. Default: \code{NULL} (single group).
#' @param ncol Integer. Number of columns in the faceted output.
#'   Default: \code{2}.
#'
#' @return A \code{ggplot} object.
#'
#' @rdname PlotQCMetrics
#' @export
setMethod("PlotQCMetrics", "MatisseObject",
          function(object, features = NULL, group_by = NULL, ncol = 2L) {
  meta <- object@isoform_metadata

  if (is.null(features)) {
    features <- colnames(meta)[vapply(meta, is.numeric, logical(1))]
  }

  missing <- setdiff(features, colnames(meta))
  if (length(missing) > 0) {
    rlang::abort(paste0(
      "Features not found in isoform_metadata: ",
      paste(missing, collapse = ", ")))
  }

  cells <- .get_cells(object)
  df    <- meta[cells, features, drop = FALSE]
  df$cell <- rownames(df)
  df$group <- if (!is.null(group_by))
    .get_seurat_meta_col(object, group_by) else "All cells"

  df_long <- tidyr::pivot_longer(
    df,
    cols      = dplyr::all_of(features),
    names_to  = "metric",
    values_to = "value"
  )

  ggplot2::ggplot(df_long, ggplot2::aes(
    x = .data$group, y = .data$value, fill = .data$group)) +
    ggplot2::geom_violin(trim = FALSE, scale = "width") +
    ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
    ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = ncol) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
    .matisse_theme()
})

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.require_psi <- function(object) {
  if (is.null(object@seurat[["psi"]])) {
    rlang::abort("PSI assay is NULL. Run CalculatePSI() first.")
  }
}

.require_event <- function(object, event_id) {
  psi_cx <- GetPSI(object)
  if (!event_id %in% colnames(psi_cx)) {
    rlang::abort(paste0("Event '", event_id,
                        "' not found in PSI assay feature names."))
  }
}

.get_seurat_meta_col <- function(object, col) {
  if (is.null(object@seurat)) rlang::abort("No Seurat object embedded.")
  meta <- object@seurat@meta.data
  if (!col %in% colnames(meta)) {
    rlang::abort(paste0("Column '", col, "' not found in Seurat meta.data."))
  }
  meta[[col]]
}

`%||%` <- function(x, y) if (!is.null(x)) x else y
