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
# PlotUMAP
# ---------------------------------------------------------------------------

#' UMAP plot coloured by any feature
#'
#' Overlays the value of a feature on the UMAP embedding stored in the
#' embedded Seurat object. Pass a PSI event ID to colour by splicing ratio,
#' a junction ID to colour by junction read counts, or a gene name to colour
#' by gene expression.
#'
#' @param object A \code{MatisseObject} with a UMAP reduction.
#' @param feature Character. Feature to plot. Can be a PSI event ID (e.g.
#'   \code{"SE:chr1:100-200:300-400:+"}), a junction ID, or a gene name.
#' @param reduction Character. Name of the dimensionality reduction to use.
#'   Default: \code{"umap"}.
#' @param dims Integer vector of length 2 selecting which dimensions to plot.
#'   Default: \code{c(1, 2)}.
#' @param pt_size Numeric. Point size. Default: \code{0.5}.
#' @param na_colour Character. Colour for cells with no data.
#'   Default: \code{"grey80"}.
#' @param title Character. Plot title. Defaults to the feature name.
#'
#' @return A \code{ggplot} object.
#'
#' @rdname PlotUMAP
#' @export
setMethod("PlotUMAP", "MatisseObject",
          function(object, feature,
                   reduction  = "umap",
                   dims       = c(1L, 2L),
                   pt_size    = 0.5,
                   na_colour  = "grey80",
                   title      = NULL) {
  vals <- .get_feature_values(object, feature)

  emb       <- SeuratObject::Embeddings(object@seurat, reduction = reduction)
  dim_names <- colnames(emb)[dims]

  df <- data.frame(
    x   = emb[, dims[1]],
    y   = emb[, dims[2]],
    val = vals
  )

  plot_title <- title %||% feature

  ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y,
                                    colour = .data$val)) +
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
# PlotViolin
# ---------------------------------------------------------------------------

#' Violin plot of feature values split by cell group
#'
#' @param object A \code{MatisseObject}.
#' @param feature Character. Feature to plot (PSI event ID, junction ID, or
#'   gene name).
#' @param group_by Character. Column in \code{Seurat::meta.data} to split
#'   cells by. Default: \code{"seurat_clusters"}.
#' @param colours Named character vector mapping group levels to colours.
#'   Default: \code{NULL} (uses ggplot2 defaults).
#' @param add_points Logical. Overlay individual cell values as jittered
#'   points. Default: \code{FALSE}.
#' @param title Character. Plot title. Defaults to the feature name.
#'
#' @return A \code{ggplot} object.
#'
#' @rdname PlotViolin
#' @export
setMethod("PlotViolin", "MatisseObject",
          function(object, feature,
                   group_by   = "seurat_clusters",
                   colours    = NULL,
                   add_points = FALSE,
                   title      = NULL) {
  vals       <- .get_feature_values(object, feature)
  group_vals <- .get_seurat_meta_col(object, group_by)

  df <- data.frame(val = vals, group = group_vals)
  df <- df[!is.na(df$val), ]

  plot_title <- title %||% feature

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data$group, y = .data$val, fill = .data$group)) +
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
# PlotHeatmap
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
#' @rdname PlotHeatmap
#' @export
setMethod("PlotHeatmap", "MatisseObject",
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
# PlotCoverage
# ---------------------------------------------------------------------------

#' Junction coverage bar plot for a gene
#'
#' Aggregates junction read counts across selected cells and plots a bar chart
#' ordered by genomic position. Only available for objects in junction mode.
#'
#' @param object A \code{MatisseObject} in junction mode with
#'   \code{junction_data} populated.
#' @param gene Character. Gene name to plot.
#' @param cells Character vector of cell barcodes to aggregate over.
#'   Default: all cells.
#' @param log_scale Logical. Use log10 y-axis. Default: \code{FALSE}.
#'
#' @return A \code{ggplot} object.
#'
#' @rdname PlotCoverage
#' @export
setMethod("PlotCoverage", "MatisseObject",
          function(object, gene, cells = NULL, log_scale = FALSE) {
  jxn_counts <- GetJunctionCounts(object)
  if (is.null(jxn_counts)) {
    rlang::abort(
      "No junction assay found. PlotCoverage() requires a junction-mode object ",
      "created with `junction_counts`.")
  }
  if (nrow(object@junction_data) == 0) {
    rlang::abort("junction_data is empty. Pass `junction_data` to CreateMatisseObject().")
  }

  jd        <- object@junction_data
  gene_jxns <- jd$junction_id[jd$gene_id == gene]
  if (length(gene_jxns) == 0) {
    rlang::abort(paste0("No junctions found for gene: ", gene))
  }

  present <- intersect(gene_jxns, colnames(jxn_counts))
  if (length(present) == 0) {
    rlang::abort(paste0("No junctions for '", gene,
                        "' are present in the junction assay."))
  }

  sub_cells <- cells %||% .get_cells(object)
  sub_mat   <- jxn_counts[sub_cells, present, drop = FALSE]
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
#' @param object A \code{MatisseObject} with QC metrics computed by
#'   \code{\link{ComputeIsoformQC}}.
#' @param features Character vector of QC metric names. Must be columns in
#'   \code{MatisseMeta(object)}. Default: all numeric columns.
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
  meta <- MatisseMeta(object)   # now returns seurat@meta.data

  if (is.null(features)) {
    features <- colnames(meta)[vapply(meta, is.numeric, logical(1))]
  }

  missing <- setdiff(features, colnames(meta))
  if (length(missing) > 0) {
    rlang::abort(paste0(
      "Features not found in metadata: ",
      paste(missing, collapse = ", "),
      ". Run ComputeIsoformQC() first."))
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
  if (is.null(.get_assay_safe(object@seurat, "psi"))) {
    rlang::abort("PSI assay is NULL. Run CalculatePSI() first.")
  }
}

# Retrieve values for a feature: checks PSI first, then junction counts,
# then gene expression from the default assay.
.get_feature_values <- function(object, feature) {
  cells <- .get_cells(object)

  # 1. PSI (most common in both modes after CalculatePSI)
  psi_cx <- GetPSI(object)
  if (!is.null(psi_cx) && feature %in% colnames(psi_cx)) {
    return(as.numeric(psi_cx[cells, feature]))
  }

  # 2. Junction counts (junction mode)
  jxn <- GetJunctionCounts(object)
  if (!is.null(jxn) && feature %in% colnames(jxn)) {
    return(as.numeric(jxn[cells, feature]))
  }

  # 3. Gene expression from default assay
  if (!is.null(object@seurat)) {
    expr <- tryCatch(
      SeuratObject::GetAssayData(object@seurat, layer = "data"),
      error = function(e) NULL
    )
    if (!is.null(expr) && feature %in% rownames(expr)) {
      return(as.numeric(expr[feature, cells]))
    }
  }

  # 4. Fall back to no PSI assay error if PSI is NULL
  if (is.null(psi_cx)) {
    rlang::abort("PSI assay is NULL. Run CalculatePSI() first.")
  }
  rlang::abort(paste0(
    "Feature '", feature, "' not found in PSI events, junction IDs, ",
    "or gene expression."))
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
