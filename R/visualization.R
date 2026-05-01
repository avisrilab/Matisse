#' @include MatisseObject-class.R
#' @include generics.R
NULL

# Shared ggplot2 theme for Matisse plots.
# Match Seurat's DimPlot / FeaturePlot look: cowplot::theme_cowplot() is the
# base Seurat uses internally (via SingleDimPlot), and Seurat::CenterTitle()
# centers the plot title the way DimPlot does. This keeps Matisse plots feeling
# native alongside Seurat plots in the same notebook.
.matisse_theme <- function() {
  cowplot::theme_cowplot() +
    Seurat::CenterTitle() +
    ggplot2::theme(
      plot.subtitle   = ggplot2::element_text(hjust = 0.5, colour = "grey40"),
      legend.position = "right"
    )
}

# ---------------------------------------------------------------------------
# PlotUMAP
# ---------------------------------------------------------------------------

#' UMAP plot -- by group (DimPlot-style) or by feature (FeaturePlot-style)
#'
#' Two modes, switched by \code{feature}:
#' \itemize{
#'   \item \strong{Group mode} (default, \code{feature = NULL}): cells coloured
#'     discretely by a metadata column -- like Seurat's \code{DimPlot}. Pass
#'     \code{label = TRUE} to add group labels at cluster centroids.
#'   \item \strong{Feature mode} (\code{feature} supplied): cells coloured on a
#'     continuous scale by a PSI event ID, junction or transcript ID, gene
#'     name, or numeric metadata column -- like Seurat's \code{FeaturePlot}.
#'     The palette adapts to the feature type: diverging RdBu \[0,1\] for
#'     PSI; sequential for counts and expression.
#' }
#'
#' @param object A \code{MatisseObject} with a UMAP reduction.
#' @param feature Character. Optional feature to plot -- a PSI event ID (e.g.
#'   \code{"SE:chr1:100-200:300-400:+"}), a junction or transcript ID, a gene
#'   name, or a numeric cell-metadata column. \code{NULL} (the default)
#'   selects group mode.
#' @param reduction Character. Name of the dimensionality reduction to use.
#'   Default: \code{"umap"}.
#' @param dims Integer vector of length 2 selecting which dimensions to plot.
#'   Default: \code{c(1, 2)}.
#' @param group_by Character. Metadata column to colour by in group mode.
#'   Default: \code{"seurat_clusters"}. Ignored in feature mode.
#' @param pt_size Numeric. Point size. Default: \code{0.5}.
#' @param na_colour Character. Colour for cells with no data.
#'   Default: \code{"grey80"}.
#' @param label Logical. Group-mode only -- add group labels at cluster
#'   centroids. Default: \code{FALSE}.
#' @param title Character. Plot title. Defaults to the feature name in
#'   feature mode and \code{NULL} in group mode.
#'
#' @return A \code{ggplot} object.
#'
#' @rdname PlotUMAP
#' @export
setMethod("PlotUMAP", "MatisseObject",
          function(object, feature = NULL,
                   reduction  = "umap",
                   dims       = c(1L, 2L),
                   group_by   = "seurat_clusters",
                   pt_size    = 0.5,
                   na_colour  = "grey80",
                   label      = FALSE,
                   title      = NULL, ...) {

  if (isTRUE(label) && !is.null(feature)) {
    rlang::abort(paste0(
      "`label = TRUE` is only valid in group mode (feature = NULL). ",
      "Centroid labels are not meaningful for a continuous colour scale."))
  }

  # ----- Group mode (DimPlot-style) -----
  # Resolve the grouping vector first so a missing metadata column errors
  # before we touch the embedding.
  if (is.null(feature)) {
    vals <- .get_seurat_meta_col(object, group_by)
    if (!is.factor(vals) && !is.character(vals)) vals <- factor(vals)

    emb       <- SeuratObject::Embeddings(object@seurat, reduction = reduction)
    dim_names <- colnames(emb)[dims]
    df <- data.frame(
      x   = emb[, dims[1]],
      y   = emb[, dims[2]],
      val = vals
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y,
                                           colour = .data$val)) +
      ggplot2::geom_point(size = pt_size, na.rm = TRUE) +
      ggplot2::labs(
        title  = title,
        x      = dim_names[1],
        y      = dim_names[2],
        colour = group_by
      ) +
      .matisse_theme()

    if (isTRUE(label)) {
      cents <- stats::aggregate(cbind(x, y) ~ val, data = df, FUN = mean)
      p <- p + ggplot2::geom_text(
        data = cents,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$val),
        colour = "black", fontface = "bold", show.legend = FALSE)
    }
    return(p)
  }

  # ----- Feature mode (FeaturePlot-style) -----
  # Classify before embedding lookup so feature-not-found / PSI-NULL surfaces
  # a feature-specific error instead of a missing-reduction one.
  classified <- .classify_feature(object, feature)
  vals       <- classified$values

  emb       <- SeuratObject::Embeddings(object@seurat, reduction = reduction)
  dim_names <- colnames(emb)[dims]

  df <- data.frame(
    x   = emb[, dims[1]],
    y   = emb[, dims[2]],
    val = vals
  )

  plot_title <- title %||% feature

  # Pick a colour scale appropriate for the feature type. PSI is bounded in
  # [0,1] and benefits from a diverging palette centred on 0.5; counts and
  # expression are non-negative and use a sequential palette without a fixed
  # upper limit.
  scale_layer <- switch(classified$type,
    psi = ggplot2::scale_colour_gradientn(
      colours  = c("#2166ac", "#f7f7f7", "#d6604d"),
      na.value = na_colour,
      limits   = c(0, 1),
      name     = "PSI"
    ),
    counts = ggplot2::scale_colour_gradientn(
      colours  = c("#fff5eb", "#fd8d3c", "#7f2704"),
      na.value = na_colour,
      name     = "counts"
    ),
    expression = ggplot2::scale_colour_gradientn(
      colours  = c("#f7fcb9", "#41ab5d", "#00441b"),
      na.value = na_colour,
      name     = "expression"
    ),
    ggplot2::scale_colour_gradientn(
      colours  = c("#f7fcb9", "#41ab5d", "#00441b"),
      na.value = na_colour,
      name     = "value"
    )
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y,
                                    colour = .data$val)) +
    ggplot2::geom_point(size = pt_size, na.rm = TRUE) +
    scale_layer +
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

#' Violin plot of one or more features split by cell group
#'
#' Plots any feature -- PSI event ID, junction ID, gene name, or metadata
#' column -- as a violin/box split by a cell-level grouping variable.
#' When \code{feature = NULL} (the default), automatically plots the standard
#' QC metrics \code{nCount_isoform}, \code{nFeature_isoform}, and
#' \code{nPercent_isoform} as a faceted panel.
#' Pass a character vector to \code{feature} to plot multiple features in
#' facets with free y-axis scales.
#'
#' @param object A \code{MatisseObject}.
#' @param feature Character scalar or vector of features to plot (PSI event
#'   IDs, junction IDs, gene names, or metadata column names). \code{NULL}
#'   (the default) selects the standard QC metrics automatically.
#' @param group_by Character. Column in \code{Seurat::meta.data} to split
#'   cells by. Default: \code{"seurat_clusters"}.
#' @param colours Named character vector mapping group levels to colours.
#'   Default: \code{NULL} (uses ggplot2 defaults).
#' @param add_points Logical. Overlay individual cell values as jittered
#'   points. Default: \code{FALSE}.
#' @param title Character. Plot title. Defaults to the feature name when a
#'   single feature is requested; \code{NULL} for multiple features.
#' @param ncol Integer. Number of facet columns when \code{feature} is a
#'   vector of length > 1. Default: \code{2}.
#' @param ... Ignored; included for S4 generic compatibility.
#'
#' @return A \code{ggplot} object.
#'
#' @rdname PlotViolin
#' @export
setMethod("PlotViolin", "MatisseObject",
          function(object, feature = NULL,
                   group_by   = "seurat_clusters",
                   colours    = NULL,
                   add_points = FALSE,
                   title      = NULL,
                   ncol       = 2L, ...) {
  # Auto-select QC metrics when no feature is specified
  if (is.null(feature)) {
    candidates <- c("nCount_isoform", "nFeature_isoform", "nPercent_isoform")
    feature    <- intersect(candidates, colnames(MatisseMeta(object)))
    if (length(feature) == 0L) {
      rlang::abort(paste0(
        "No QC metrics found in cell metadata. ",
        "Run CreateMatisseObject() and CalculatePSI() first."))
    }
  }

  # Resolve feature values first -- gives the right error when PSI is NULL
  vals_list  <- lapply(feature, function(f) .get_feature_values(object, f))
  group_vals <- .get_seurat_meta_col(object, group_by)

  # Build long-form data frame -- one row per (cell, feature)
  df <- do.call(rbind, mapply(function(f, vals) {
    data.frame(val = vals, group = group_vals, feature_name = f,
               stringsAsFactors = FALSE)
  }, feature, vals_list, SIMPLIFY = FALSE))
  df <- df[!is.na(df$val), ]

  plot_title <- title %||% if (length(feature) == 1L) feature else NULL

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data$group, y = .data$val, fill = .data$group)) +
    ggplot2::geom_violin(trim = FALSE, scale = "width") +
    ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
    ggplot2::labs(title = plot_title, x = group_by, y = NULL) +
    ggplot2::theme(legend.position = "none") +
    .matisse_theme()

  if (length(feature) > 1L) {
    p <- p + ggplot2::facet_wrap(~ feature_name, scales = "free_y", ncol = ncol)
  }

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

#' Heatmap of PSI values (events x cells)
#'
#' Draws a DoHeatmap-style tile plot with splice events on the y-axis and cells
#' on the x-axis. Events are clustered by PSI profile. When \code{group_by} is
#' supplied, cells are ordered by group and labelled with facet strips.
#'
#' @param object A \code{MatisseObject} with a \code{"psi"} assay.
#' @param events Character vector of event IDs to include.
#'   Default: \code{NULL} (top-variance events up to \code{max_events}).
#' @param cells Character vector of cell barcodes to include.
#'   Default: \code{NULL} (random sample up to \code{max_cells}).
#' @param group_by Character. Column in \code{Seurat::meta.data} used to
#'   order and label cells. Default: \code{NULL}.
#' @param max_cells Integer. Downsample to this many cells before plotting.
#'   Default: \code{500}.
#' @param max_events Integer. Cap on events to plot. When the candidate set
#'   exceeds this, the top-variance events are selected automatically.
#'   Default: \code{200}.
#' @param na_colour Character. Colour for \code{NA} entries.
#'   Default: \code{"grey90"}.
#' @param title Character. Plot title. Default: \code{"PSI Heatmap"}.
#' @param ... Ignored; included for S4 generic compatibility.
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
                   na_colour  = "grey90",
                   title      = NULL, ...) {
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

  # Cap events by variance before densifying
  if (length(events) > max_events) {
    psi_for_var <- psi_cx[cells, events, drop = FALSE]
    col_vars    <- apply(psi_for_var, 2, stats::var, na.rm = TRUE)
    col_vars[is.na(col_vars)] <- 0
    events <- events[order(col_vars, decreasing = TRUE)[seq_len(max_events)]]
    rlang::inform(paste0(
      "Showing top ", max_events, " highest-variance events. ",
      "Pass `events` explicitly or increase `max_events` to change this."))
  }

  # Densify: at most max_cells x max_events
  psi_sub <- as.matrix(psi_cx[cells, events, drop = FALSE])

  # Cluster events by PSI profile
  finite_mask <- is.finite(psi_sub)
  if (sum(finite_mask) > 0 && length(events) > 1L) {
    psi_for_clust               <- psi_sub
    psi_for_clust[!finite_mask] <- 0.5
    event_order <- tryCatch({
      hc <- stats::hclust(stats::dist(t(psi_for_clust)))
      hc$order
    }, error = function(e) seq_len(ncol(psi_sub)))
  } else {
    event_order <- seq_len(ncol(psi_sub))
  }

  # Order cells by group if requested
  grp <- NULL
  if (!is.null(group_by)) {
    grp   <- .get_seurat_meta_col(object, group_by)[match(cells, all_cells)]
    ord   <- order(grp)
    cells <- cells[ord]
    grp   <- grp[ord]
  }

  # Build tidy data frame: events on y, cells on x
  df <- data.frame(
    cell  = rep(cells,               times = length(events)),
    event = rep(events[event_order], each  = length(cells)),
    psi   = as.vector(psi_sub[cells, events[event_order]])
  )
  df$cell  <- factor(df$cell,  levels = cells)
  # Reverse so the first clustered event sits at the top of the y-axis
  df$event <- factor(df$event, levels = rev(events[event_order]))

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x    = .data$cell,
    y    = .data$event,
    fill = .data$psi)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colours  = c("#2166ac", "#f7f7f7", "#d6604d"),
      na.value = na_colour,
      limits   = c(0, 1),
      name     = "PSI"
    ) +
    ggplot2::labs(title = title %||% "PSI Heatmap") +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      # Adapt y-axis label size to event count: tiny size makes labels
      # illegible at high event counts but is fine when there are few rows.
      axis.text.y  = ggplot2::element_text(
        size = max(3, min(10, round(180 / max(length(events), 1L))))
      ),
      axis.title.y = ggplot2::element_blank()
    ) +
    .matisse_theme()

  # Group annotation strips (DoHeatmap style)
  if (!is.null(group_by)) {
    df$group <- factor(rep(grp, times = length(events)), levels = unique(grp))
    p <- p + ggplot2::facet_grid(. ~ group, scales = "free_x", space = "free") +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = "grey85", colour = NA),
        strip.text       = ggplot2::element_text(size = 8L, face = "bold"),
        panel.spacing    = ggplot2::unit(0.5, "mm")
      )
  }

  p
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
#' with its own read count. In \strong{transcript mode} the SE event_id is
#' parsed to derive junction coordinates; inclusion and exclusion counts come
#' from the \code{"counts"} and \code{"exclusion"} layers of the PSI assay.
#'
#' Supported event types in transcript mode: \strong{SE} (skipped exon) and
#' \strong{RI} (retained intron). Junction mode supports all event types
#' since coordinates are auto-parsed from junction IDs.
#'
#' \strong{Note on transcript-mode SE arcs:} transcript-level counting
#' aggregates reads to events, not to individual junctions. The total
#' inclusion-event count is therefore split evenly across the two SE
#' inclusion arcs in the plot. The two arcs may have differed in reality;
#' use junction-mode input if you need per-junction read counts.
#'
#' @param object A \code{MatisseObject} with a PSI assay computed.
#' @param event_id Character. A single event ID. Run
#'   \code{rownames(GetSeurat(obj)[["psi"]])} to list available IDs;
#'   typical SE format is \code{"SE:chr1:1201-2999:3201-4999:+"}.
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

  # Read event annotation from the PSI assay's meta.features (post-P1).
  # Pre-CalculatePSI, fall back to the @misc staging area.
  psi_assay <- .get_assay_safe(object@seurat, "psi")
  ed <- if (!is.null(psi_assay)) {
    mf <- psi_assay[[]]
    if (!is.null(mf) && ncol(mf) > 0L) {
      mf$event_id <- rownames(mf)
      mf
    } else NULL
  } else object@misc[["event_data"]]
  if (is.null(ed) || !event_id %in% ed$event_id) {
    rlang::abort(paste0("'", event_id, "' not found in event annotation."))
  }
  ev        <- ed[ed$event_id == event_id, , drop = FALSE]
  all_cells <- .get_cells(object)
  sub_cells <- cells %||% all_cells

  # Split cells by group; drop empty groups so .cov_counts is never asked
  # to sum over zero rows (which would produce NaN arc heights).
  if (!is.null(group_by)) {
    grp_vec        <- .get_seurat_meta_col(object, group_by)
    names(grp_vec) <- all_cells
    groups         <- split(sub_cells, grp_vec[sub_cells])
    groups         <- groups[lengths(groups) > 0L]
    if (length(groups) == 0L) {
      rlang::abort("No cells in any group of `group_by`.")
    }
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
  if (object@input.mode == "junction") {
    inc <- strsplit(ev$inclusion_features, ";", fixed = TRUE)[[1]]
    exc <- strsplit(ev$exclusion_features, ";", fixed = TRUE)[[1]]
    ids <- c(inc, exc)
    # Per-junction coordinates live in the isoform assay's meta.features,
    # auto-derived from junction IDs at construction.
    iso <- .get_assay_safe(object@seurat, "isoform")
    if (is.null(iso)) {
      rlang::abort(
        "No 'isoform' assay found; cannot draw junction coordinates.")
    }
    jd  <- iso[[]]
    idx <- match(ids, rownames(jd))
    if (anyNA(idx)) {
      missing <- ids[is.na(idx)]
      rlang::abort(paste0(
        "Junctions referenced by event '", ev$event_id,
        "' not found in junction count matrix: ",
        paste(missing, collapse = ", ")))
    }
    bad_coords <- is.na(jd$chr[idx])
    if (any(bad_coords)) {
      rlang::abort(paste0(
        "Could not parse coordinates from these junction IDs (needed for ",
        "sashimi plot): ", paste(ids[bad_coords], collapse = ", "),
        ". Supported formats: chr-start-end-strand, chr:start-end:strand, ",
        "chr_start_end_strand."))
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
        "PlotSashimi() in transcript mode does not yet support '", etype,
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
  if (object@input.mode == "junction") {
    iso <- .get_assay_safe(object@seurat, "isoform")
    jxn <- if (!is.null(iso))
      Matrix::t(.get_assay_layer(iso, "counts")) else NULL
    inc <- strsplit(ev$inclusion_features, ";", fixed = TRUE)[[1]]
    exc <- strsplit(ev$exclusion_features, ";", fixed = TRUE)[[1]]
    ids <- intersect(c(inc, exc), colnames(jxn))
    tot <- as.numeric(Matrix::colSums(jxn[cells, ids, drop = FALSE]))
    data.frame(junction_id = ids, count = tot, stringsAsFactors = FALSE)
  } else {
    .require_psi(object)
    eid       <- ev$event_id
    psi_assay <- .get_assay_safe(object@seurat, "psi")
    inc_cx    <- if (!is.null(psi_assay))
      Matrix::t(.get_assay_layer(psi_assay, "counts")) else NULL
    exc_layer <- if (!is.null(psi_assay))
      .get_assay_layer(psi_assay, "exclusion") else NULL
    exc_cx    <- if (!is.null(exc_layer) && length(exc_layer) > 0)
      Matrix::t(exc_layer) else NULL
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

# Map an arc_scale label to the corresponding R function used to scale
# arc heights by read count. Single source of truth so .cov_arc_paths and
# .cov_draw stay aligned if a new scale option is ever added.
.arc_scale_fn <- function(arc_scale) {
  switch(arc_scale,
    sqrt   = sqrt,
    linear = identity,
    log    = log1p
  )
}

# Build arc path points (n_pts per arc) and return long data.frame
.cov_arc_paths <- function(arc_data, arc_scale) {
  scale_fn <- .arc_scale_fn(arc_scale)
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
  scale_fn <- .arc_scale_fn(arc_scale)

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

  if (facet) p <- p + ggplot2::facet_wrap(~ group, ncol = 1L)
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

# Classify a feature by source — used by PlotUMAP to pick the right colour
# scale. Returns list(values, type) where type is one of
# "psi" / "counts" / "expression" / "metadata".
.classify_feature <- function(object, feature) {
  cells <- .get_cells(object)

  # 1. PSI event
  psi_cx <- GetPSI(object)
  if (!is.null(psi_cx) && feature %in% colnames(psi_cx)) {
    return(list(values = as.numeric(psi_cx[cells, feature]), type = "psi"))
  }

  # 2. Junction counts (junction mode)
  iso <- if (object@input.mode == "junction")
    .get_assay_safe(object@seurat, "isoform") else NULL
  jxn <- if (!is.null(iso))
    Matrix::t(.get_assay_layer(iso, "counts")) else NULL
  if (!is.null(jxn) && feature %in% colnames(jxn)) {
    return(list(values = as.numeric(jxn[cells, feature]), type = "counts"))
  }

  # 3. Gene expression from default assay
  if (!is.null(object@seurat)) {
    expr <- tryCatch(
      SeuratObject::GetAssayData(object@seurat, layer = "data"),
      error = function(e) NULL
    )
    if (!is.null(expr) && feature %in% rownames(expr)) {
      return(list(values = as.numeric(expr[feature, cells]),
                  type   = "expression"))
    }
  }

  # 4. Seurat cell metadata
  if (!is.null(object@seurat) &&
      feature %in% colnames(object@seurat@meta.data)) {
    return(list(values = as.numeric(object@seurat@meta.data[cells, feature]),
                type   = "metadata"))
  }

  # Not found — fall through to .get_feature_values for the abort message.
  list(values = .get_feature_values(object, feature), type = "metadata")
}

# Retrieve values for a feature: checks PSI first, then junction counts,
# then gene expression from the default assay. Used by PlotViolin and
# friends; PlotUMAP uses .classify_feature so it can pick a colour scale
# matched to the feature type.
.get_feature_values <- function(object, feature) {
  cells  <- .get_cells(object)
  psi_cx <- GetPSI(object)

  # 1. PSI (most common in both modes after CalculatePSI)
  if (!is.null(psi_cx) && feature %in% colnames(psi_cx)) {
    return(as.numeric(psi_cx[cells, feature]))
  }

  # 2. Junction counts (junction mode) — read straight from the isoform assay
  iso <- if (object@input.mode == "junction")
    .get_assay_safe(object@seurat, "isoform") else NULL
  jxn <- if (!is.null(iso))
    Matrix::t(.get_assay_layer(iso, "counts")) else NULL
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

  # 4. Seurat metadata (nCount_isoform, nPercent_isoform, cell_type, etc.)
  if (!is.null(object@seurat) &&
      feature %in% colnames(object@seurat@meta.data)) {
    return(as.numeric(object@seurat@meta.data[cells, feature]))
  }

  # 5. Not found
  if (is.null(psi_cx)) {
    rlang::abort("PSI assay is NULL. Run CalculatePSI() first.")
  }
  rlang::abort(paste0(
    "Feature '", feature, "' not found in PSI events, junction IDs, ",
    "gene expression, or cell metadata."))
}

.get_seurat_meta_col <- function(object, col) {
  if (is.null(object@seurat)) rlang::abort("No Seurat object embedded.")
  meta <- object@seurat@meta.data
  if (!col %in% colnames(meta)) {
    rlang::abort(paste0("Column '", col, "' not found in Seurat meta.data."))
  }
  meta[[col]]
}
