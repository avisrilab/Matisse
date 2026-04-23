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
                   title      = NULL, ...) {
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
                   title      = NULL, ...) {
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
#'   Default: \code{NULL} (top-variance events up to \code{max_events}).
#' @param cells Character vector of cell barcodes to include.
#'   Default: \code{NULL} (random sample up to \code{max_cells}).
#' @param group_by Character. Column in \code{Seurat::meta.data} used to
#'   annotate and order cells. Default: \code{NULL}.
#' @param max_cells Integer. Downsample to this many cells before plotting.
#'   Default: \code{500}.
#' @param max_events Integer. Cap on events to plot. When the candidate set
#'   exceeds this, the top-variance events are selected automatically.
#'   Default: \code{200}.
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
                   max_events = 200L,
                   na_colour  = "grey90", ...) {
  .require_psi(object)

  psi_cx     <- GetPSI(object)  # cells x events (sparse)
  all_cells  <- .get_cells(object)
  all_events <- colnames(psi_cx)

  # Validate and intersect user-supplied vectors
  cells  <- cells  %||% all_cells
  events <- events %||% all_events

  bad_cells  <- setdiff(cells,  all_cells)
  bad_events <- setdiff(events, all_events)
  if (length(bad_cells)  > 0) rlang::warn("Some cell barcodes not found; skipping.")
  if (length(bad_events) > 0) rlang::warn("Some event IDs not found; skipping.")
  cells  <- intersect(cells,  all_cells)
  events <- intersect(events, all_events)

  # Cap cells first (cheap on sparse matrix)
  if (length(cells) > max_cells) cells <- sample(cells, max_cells)

  # Cap events by variance — must happen BEFORE as.matrix() to avoid huge allocs
  if (length(events) > max_events) {
    # Compute per-event variance on the sparse submatrix (still sparse here)
    psi_for_var <- psi_cx[cells, events, drop = FALSE]
    col_vars    <- apply(psi_for_var, 2, stats::var, na.rm = TRUE)
    col_vars[is.na(col_vars)] <- 0
    events <- events[order(col_vars, decreasing = TRUE)[seq_len(max_events)]]
    rlang::inform(paste0(
      "Showing top ", max_events, " highest-variance events. ",
      "Pass `events` explicitly or increase `max_events` to change this."))
  }

  # Now safe to densify: at most max_cells x max_events
  psi_sub <- as.matrix(psi_cx[cells, events, drop = FALSE])

  # Cluster events by PSI profile (small matrix now — safe)
  finite_mask <- is.finite(psi_sub)
  if (sum(finite_mask) > 0 && length(events) > 1L) {
    psi_for_clust              <- psi_sub
    psi_for_clust[!finite_mask] <- 0.5
    event_order <- tryCatch({
      hc <- stats::hclust(stats::dist(t(psi_for_clust)))
      hc$order
    }, error = function(e) seq_len(ncol(psi_sub)))
  } else {
    event_order <- seq_len(ncol(psi_sub))
  }

  # Order cells by group if requested
  if (!is.null(group_by)) {
    grp   <- .get_seurat_meta_col(object, group_by)[match(cells, all_cells)]
    cells <- cells[order(grp)]
  }

  df <- data.frame(
    cell  = rep(cells,               each  = length(events)),
    event = rep(events[event_order], times = length(cells)),
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
          function(object, gene, cells = NULL, log_scale = FALSE, ...) {
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
          function(object, features = NULL, group_by = NULL, ncol = 2L, ...) {
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
# PlotSashimi
# ---------------------------------------------------------------------------

#' Sashimi-style coverage plot for a splice event
#'
#' Draws junction arcs scaled by aggregate read count over a schematic gene
#' structure. Arcs are coloured by role: inclusion (blue) vs exclusion (red).
#'
#' In \strong{junction mode} each arc corresponds to an individual junction
#' with its own read count. In \strong{event mode} the SE event_id is parsed
#' to derive junction coordinates; inclusion and exclusion counts come from the
#' \code{"counts"} and \code{"exclusion"} layers of the PSI assay.
#'
#' Supported event types in event mode: \strong{SE} (skipped exon) and
#' \strong{RI} (retained intron). Junction mode supports all event types
#' since coordinates come directly from \code{junction_data}.
#'
#' @param object A \code{MatisseObject} with a PSI assay computed.
#' @param event_id Character. Event ID as stored in \code{event_data}, e.g.
#'   \code{"SE:chr1:1201-2999:3201-4999:+"}.
#' @param cells Character vector of cell barcodes to aggregate over.
#'   Default: all cells.
#' @param group_by Character. Column in Seurat meta.data to facet by.
#'   Default: \code{NULL} (all cells pooled).
#' @param arc_scale Character. How to scale arc height to read count:
#'   \code{"sqrt"} (default), \code{"linear"}, or \code{"log"}.
#' @param colours Named character vector with elements \code{"inclusion"} and
#'   \code{"exclusion"} giving arc colours.
#' @param title Character. Plot title. Defaults to \code{event_id}.
#'
#' @return A \code{ggplot} object.
#'
#' @rdname PlotSashimi
#' @export
setMethod("PlotSashimi", "MatisseObject",
          function(object, event_id,
                   cells     = NULL,
                   group_by  = NULL,
                   arc_scale = c("sqrt", "linear", "log"),
                   colours   = c(inclusion = "#4393c3", exclusion = "#d6604d"),
                   title     = NULL, ...) {
  arc_scale <- match.arg(arc_scale)

  ed <- object@event_data
  if (!event_id %in% ed$event_id) {
    rlang::abort(paste0("'", event_id, "' not found in event_data."))
  }
  ev        <- ed[ed$event_id == event_id, , drop = FALSE]
  all_cells <- .get_cells(object)
  sub_cells <- cells %||% all_cells

  # Split cells by group
  if (!is.null(group_by)) {
    grp_vec        <- .get_seurat_meta_col(object, group_by)
    names(grp_vec) <- all_cells
    groups         <- split(sub_cells, grp_vec[sub_cells])
  } else {
    groups <- list(All = sub_cells)
  }

  jxn_coords <- .cov_jxn_coords(object, ev)

  arc_data <- do.call(rbind, lapply(names(groups), function(g) {
    cnts     <- .cov_counts(object, ev, groups[[g]])
    df       <- merge(jxn_coords, cnts, by = "junction_id", all.x = TRUE)
    df$count <- ifelse(is.na(df$count), 0, df$count)
    df$group <- g
    df
  }))

  .cov_draw(arc_data, jxn_coords, arc_scale, colours,
            title %||% event_id, !is.null(group_by))
})

# Return data.frame: junction_id | start | end | role
.cov_jxn_coords <- function(object, ev) {
  if (object@mode == "junction") {
    inc <- strsplit(ev$inclusion_junctions, ";", fixed = TRUE)[[1]]
    exc <- strsplit(ev$exclusion_junctions, ";", fixed = TRUE)[[1]]
    ids <- c(inc, exc)
    jd  <- object@junction_data
    idx <- match(ids, jd$junction_id)
    if (anyNA(idx)) {
      missing <- ids[is.na(idx)]
      rlang::abort(paste0(
        "Junctions not found in junction_data: ",
        paste(missing, collapse = ", ")))
    }
    data.frame(
      junction_id = ids,
      chr         = jd$chr[idx],
      start       = jd$start[idx],
      end         = jd$end[idx],
      strand      = jd$strand[idx],
      role        = c(rep("inclusion", length(inc)),
                      rep("exclusion", length(exc))),
      stringsAsFactors = FALSE
    )
  } else {
    etype <- strsplit(ev$event_id, ":", fixed = TRUE)[[1L]][1L]
    switch(etype,
      SE = .cov_parse_se(ev),
      RI = .cov_parse_ri(ev),
      rlang::abort(paste0(
        "PlotSashimi() in event mode does not yet support '", etype,
        "' events. Supported types: SE, RI."))
    )
  }
}

# Parse SE event row into junction coord table (chr + strand from ev)
# Format: SE:chr:donor1-acceptor1:donor2-acceptor2:strand
.cov_parse_se <- function(ev) {
  event_id <- ev$event_id
  parts    <- strsplit(event_id, ":", fixed = TRUE)[[1L]]
  if (length(parts) < 5L || parts[1L] != "SE") {
    rlang::abort(paste0(
      "Malformed SE event ID: '", event_id, "'. ",
      "Expected format: 'SE:chr:start1-end1:start2-end2:strand'."))
  }
  chr    <- ev$chr
  strand <- ev$strand
  p3     <- as.integer(strsplit(parts[3L], "-", fixed = TRUE)[[1L]])
  p4     <- as.integer(strsplit(parts[4L], "-", fixed = TRUE)[[1L]])
  data.frame(
    junction_id = c("inc_jxn1", "inc_jxn2", "exc_jxn"),
    chr         = chr,
    start       = c(p3[1L], p4[1L], p3[1L]),
    end         = c(p3[2L], p4[2L], p4[2L]),
    strand      = strand,
    role        = c("inclusion", "inclusion", "exclusion"),
    stringsAsFactors = FALSE
  )
}

# Parse RI event row into junction coord table (chr + strand from ev)
# Format: RI:chr:exon1_end:intron_start-intron_end:exon2_start:strand
# inc_jxn represents the retained intron body; exc_jxn is the normal splice.
.cov_parse_ri <- function(ev) {
  event_id <- ev$event_id
  parts    <- strsplit(event_id, ":", fixed = TRUE)[[1L]]
  if (length(parts) < 6L || parts[1L] != "RI") {
    rlang::abort(paste0(
      "Malformed RI event ID: '", event_id, "'. ",
      "Expected format: 'RI:chr:exon1_end:intron_start-intron_end:exon2_start:strand'."))
  }
  chr        <- ev$chr
  strand     <- ev$strand
  exon1_end  <- as.integer(parts[3L])
  intron     <- as.integer(strsplit(parts[4L], "-", fixed = TRUE)[[1L]])
  exon2_start <- as.integer(parts[5L])
  data.frame(
    junction_id = c("inc_jxn", "exc_jxn"),
    chr         = chr,
    start       = c(intron[1L], exon1_end),
    end         = c(intron[2L], exon2_start),
    strand      = strand,
    role        = c("inclusion", "exclusion"),
    stringsAsFactors = FALSE
  )
}

# Return data.frame: junction_id | count
.cov_counts <- function(object, ev, cells) {
  if (object@mode == "junction") {
    jxn <- GetJunctionCounts(object)
    inc <- strsplit(ev$inclusion_junctions, ";", fixed = TRUE)[[1]]
    exc <- strsplit(ev$exclusion_junctions, ";", fixed = TRUE)[[1]]
    ids <- intersect(c(inc, exc), colnames(jxn))
    tot <- as.numeric(Matrix::colSums(jxn[cells, ids, drop = FALSE]))
    data.frame(junction_id = ids, count = tot, stringsAsFactors = FALSE)
  } else {
    .require_psi(object)
    eid     <- ev$event_id
    inc_cx  <- GetInclusionCounts(object)
    exc_cx  <- GetExclusionCounts(object)
    inc_tot <- if (!is.null(inc_cx) && eid %in% colnames(inc_cx))
      sum(inc_cx[cells, eid]) else 0
    exc_tot <- if (!is.null(exc_cx) && eid %in% colnames(exc_cx))
      sum(exc_cx[cells, eid]) else 0
    etype <- strsplit(eid, ":", fixed = TRUE)[[1L]][1L]
    if (etype == "RI") {
      data.frame(
        junction_id = c("inc_jxn", "exc_jxn"),
        count       = c(inc_tot,    exc_tot),
        stringsAsFactors = FALSE
      )
    } else {
      # SE (and other 2-junction events): split inclusion evenly across the two arcs
      data.frame(
        junction_id = c("inc_jxn1", "inc_jxn2", "exc_jxn"),
        count       = c(inc_tot / 2, inc_tot / 2, exc_tot),
        stringsAsFactors = FALSE
      )
    }
  }
}

# Build arc path points (n_pts per arc) and return long data.frame
.cov_arc_paths <- function(arc_data, arc_scale) {
  scale_fn <- switch(arc_scale,
    sqrt   = sqrt,
    linear = identity,
    log    = log1p
  )
  n_pts <- 60L
  t_seq <- seq(0, pi, length.out = n_pts)
  do.call(rbind, lapply(seq_len(nrow(arc_data)), function(i) {
    row <- arc_data[i, ]
    h   <- scale_fn(max(row$count, 0))
    x1  <- row$start
    x2  <- row$end
    data.frame(
      x     = (x1 + x2) / 2 - (x2 - x1) / 2 * cos(t_seq),
      y     = h * sin(t_seq),
      role  = row$role,
      group = row$group,
      arc   = i,
      stringsAsFactors = FALSE
    )
  }))
}

# Derive exon blocks and intron backbone from junction coordinate table
.cov_gene_model <- function(jxn_coords, x_pad = 300L) {
  donors    <- sort(unique(jxn_coords$start))
  acceptors <- sort(unique(jxn_coords$end))
  x_min     <- min(donors)    - x_pad
  x_max     <- max(acceptors) + x_pad
  exon_xmin <- c(x_min,         acceptors + 1L)
  exon_xmax <- c(donors  - 1L,  x_max)
  list(
    exons   = data.frame(xmin = exon_xmin, xmax = exon_xmax,
                         ymin = -0.08, ymax = 0.08),
    intron  = data.frame(x = c(x_min, x_max), y = c(0, 0)),
    x_range = c(x_min, x_max)
  )
}

.cov_draw <- function(arc_data, jxn_coords, arc_scale, colours, title, facet) {
  scale_fn <- switch(arc_scale, sqrt = sqrt, linear = identity, log = log1p)

  arc_paths <- .cov_arc_paths(arc_data, arc_scale)
  gene      <- .cov_gene_model(jxn_coords)

  # Count labels: peak of each arc
  label_df <- arc_data
  label_df$lx <- (arc_data$start + arc_data$end) / 2
  label_df$ly <- scale_fn(pmax(arc_data$count, 0)) * 1.08

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = gene$intron,
      ggplot2::aes(x = .data$x, y = .data$y),
      colour = "grey50", linewidth = 0.4) +
    ggplot2::geom_rect(
      data = gene$exons,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                   ymin = .data$ymin, ymax = .data$ymax),
      fill = "grey55", colour = NA) +
    ggplot2::geom_path(
      data = arc_paths,
      ggplot2::aes(x      = .data$x,
                   y      = .data$y,
                   colour = .data$role,
                   group  = .data$arc),
      linewidth = 1.1, lineend = "round") +
    ggplot2::geom_text(
      data = label_df[label_df$count > 0, ],
      ggplot2::aes(x     = .data$lx,
                   y     = .data$ly,
                   label = round(.data$count),
                   colour = .data$role),
      size = 3, vjust = 0, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = colours, name = "Junction role") +
    ggplot2::scale_x_continuous(
      limits = gene$x_range,
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    ggplot2::labs(
      title = title,
      x     = paste0(jxn_coords$chr[1L], "  (", jxn_coords$strand[1L], ")"),
      y     = paste0("Reads (", arc_scale, " scaled)")
    ) +
    .matisse_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank())

  if (facet) p <- p + ggplot2::facet_wrap(~ group)
  p
}

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
