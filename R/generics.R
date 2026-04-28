#' @include MatisseObject-class.R
NULL

# ---------------------------------------------------------------------------
# Accessor generics
# ---------------------------------------------------------------------------

#' Get the embedded Seurat object
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments (unused).
#' @return A \code{Seurat} object.
#' @export
setGeneric("GetSeurat", function(object, ...) standardGeneric("GetSeurat"))

#' Get the PSI matrix
#'
#' Retrieves the PSI (Percent Spliced In) matrix from the \code{"psi"}
#' \code{Assay5} stored inside the embedded Seurat object.
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments (unused).
#' @return A sparse matrix (cells x events) of PSI values in \eqn{[0,1]}.
#'   \code{NULL} if no \code{"psi"} assay exists yet.
#' @export
setGeneric("GetPSI", function(object, ...) standardGeneric("GetPSI"))

#' Set the PSI matrix
#'
#' Replaces the \code{"data"} layer of the \code{"psi"} \code{Assay5}
#' inside the embedded Seurat object.
#'
#' @param object A \code{MatisseObject}.
#' @param value A sparse matrix (cells x events) of PSI values.
#' @return The updated \code{MatisseObject}.
#' @export
setGeneric("SetPSI", function(object, value) standardGeneric("SetPSI"))

#' Get or set cell-level metadata
#'
#' Returns the full \code{meta.data} of the embedded Seurat object, which
#' includes all per-cell QC metrics and annotations added by Matisse (e.g.
#' \code{nCount_isoform}, \code{nFeature_isoform}, \code{nPercent_isoform})
#' alongside standard Seurat columns. To add new columns, use Seurat's
#' \code{\link[SeuratObject]{AddMetaData}} (an S3 method dispatches on
#' \code{MatisseObject}).
#'
#' @param object A \code{MatisseObject}.
#' @param value A \code{data.frame} whose columns are added to cell metadata
#'   (for the setter). Rownames must match cell barcodes.
#' @param ... Additional arguments (unused).
#' @return For the getter: a \code{data.frame}. For the setter: the updated
#'   \code{MatisseObject}.
#' @export
setGeneric("MatisseMeta",
           function(object, ...) standardGeneric("MatisseMeta"))

#' @rdname MatisseMeta
#' @export
setGeneric("MatisseMeta<-",
           function(object, value) standardGeneric("MatisseMeta<-"))

# ---------------------------------------------------------------------------
# Analysis generics
# ---------------------------------------------------------------------------

#' Calculate PSI matrix from junction or transcript counts
#'
#' Typically called automatically by \code{\link{CreateMatisseObject}}; call
#' directly only to recompute with different parameters (e.g. a different
#' \code{min_coverage}).
#'
#' @param object A \code{MatisseObject}, or a sparse count matrix
#'   (cells x junctions) for raw-matrix usage.
#' @param events A \code{data.frame} defining splice events with columns
#'   \code{event_id}, \code{inclusion_features}, \code{exclusion_features}
#'   (see the method documentation for the full list of supported columns).
#' @param ... Additional arguments passed to the method.
#' @return A \code{MatisseObject} (when given one) or a PSI matrix.
#' @export
setGeneric("CalculatePSI",
           function(object, events = NULL, min_coverage = 5L,
                    verbose = TRUE, ...)
             standardGeneric("CalculatePSI"))

#' Filter cells by QC thresholds
#'
#' @param object A \code{MatisseObject}.
#' @param min_features_isoform Integer. Minimum \code{nFeature_isoform}.
#'   Default: \code{NULL}.
#' @param max_features_isoform Integer. Maximum \code{nFeature_isoform}.
#'   Default: \code{NULL}.
#' @param min_counts_isoform Integer. Minimum \code{nCount_isoform}.
#'   Default: \code{NULL}.
#' @param max_counts_isoform Integer. Maximum \code{nCount_isoform}.
#'   Default: \code{NULL}.
#' @param min_pct_isoform Numeric (0-100). Minimum \code{nPercent_isoform}.
#'   Default: \code{NULL}.
#' @param custom_filters Named list of \code{c(min, max)} bounds for arbitrary
#'   metadata columns. Use \code{NA} for one-sided bounds.
#' @param verbose Logical. Default: \code{TRUE}.
#' @return The filtered \code{MatisseObject}.
#' @export
setGeneric("FilterCells",
           function(object,
                    min_features_isoform = NULL,
                    max_features_isoform = NULL,
                    min_counts_isoform   = NULL,
                    max_counts_isoform   = NULL,
                    min_pct_isoform      = NULL,
                    custom_filters       = NULL,
                    verbose              = TRUE, ...)
             standardGeneric("FilterCells"))

#' Filter splice events by coverage or variance
#'
#' @param object A \code{MatisseObject}.
#' @param min_cells_covered Integer. Minimum cells with non-NA PSI.
#'   Default: \code{10}.
#' @param min_psi_variance Numeric. Minimum PSI variance across covered cells.
#'   Default: \code{NULL}.
#' @param verbose Logical. Default: \code{TRUE}.
#' @return The filtered \code{MatisseObject}.
#' @export
setGeneric("FilterEvents",
           function(object,
                    min_cells_covered = 10L,
                    min_psi_variance  = NULL,
                    verbose           = TRUE, ...)
             standardGeneric("FilterEvents"))

# ---------------------------------------------------------------------------
# Visualization generics
# ---------------------------------------------------------------------------

#' UMAP plot coloured by any feature
#'
#' Overlays the value of a feature on the UMAP embedding stored in the
#' embedded Seurat object. Pass a PSI event ID, a junction ID, a gene
#' name, or a cell-metadata column. The colour scale adapts automatically:
#' diverging RdBu \[0,1\] for PSI; sequential viridis for counts and
#' expression.
#'
#' @param object A \code{MatisseObject}.
#' @param feature Character. Feature to visualise.
#' @param ... Additional arguments (see \code{\link{PlotUMAP}}).
#' @return A \code{ggplot} object.
#' @export
setGeneric("PlotUMAP",
           function(object, feature,
                    reduction = "umap",
                    dims      = c(1L, 2L),
                    pt_size   = 0.5,
                    na_colour = "grey80",
                    title     = NULL, ...)
             standardGeneric("PlotUMAP"))

#' Violin plot of feature values split by cell group
#'
#' When \code{feature} is \code{NULL} (the default), automatically plots the
#' QC metrics \code{nCount_isoform}, \code{nFeature_isoform}, and
#' \code{nPercent_isoform} as a faceted panel.
#'
#' @param object A \code{MatisseObject}.
#' @param feature Character vector of features to plot (PSI event IDs,
#'   junction IDs, gene names, or metadata columns). \code{NULL} defaults to
#'   the standard QC metrics.
#' @param group_by Character. Metadata column to split cells by.
#'   Default: \code{"seurat_clusters"}.
#' @param colours Named character vector mapping group levels to colours.
#'   Default: \code{NULL}.
#' @param add_points Logical. Overlay jittered cell values. Default:
#'   \code{FALSE}.
#' @param title Character. Plot title. Default: feature name (single feature)
#'   or \code{NULL} (multiple features).
#' @param ncol Integer. Number of facet columns when \code{feature} is a
#'   vector. Default: \code{2}.
#' @return A \code{ggplot} object.
#' @export
setGeneric("PlotViolin",
           function(object, feature = NULL,
                    group_by   = "seurat_clusters",
                    colours    = NULL,
                    add_points = FALSE,
                    title      = NULL,
                    ncol       = 2L, ...)
             standardGeneric("PlotViolin"))

#' Heatmap of PSI values (events x cells, DoHeatmap style)
#'
#' @param object A \code{MatisseObject}.
#' @param events Character vector of event IDs. Default: top-variance events
#'   up to \code{max_events}.
#' @param cells Character vector of cell barcodes. Default: random sample up
#'   to \code{max_cells}.
#' @param group_by Character. Metadata column to annotate and order cells.
#'   Default: \code{NULL}.
#' @param max_cells Integer. Cell downsample cap. Default: \code{500}.
#' @param max_events Integer. Event cap; top-variance events selected when
#'   exceeded. Default: \code{200}.
#' @param na_colour Character. Colour for \code{NA} entries.
#'   Default: \code{"grey90"}.
#' @param title Character. Plot title. Default: \code{"PSI Heatmap"}.
#' @return A \code{ggplot} object.
#' @export
setGeneric("PlotHeatmap",
           function(object,
                    events     = NULL,
                    cells      = NULL,
                    group_by   = NULL,
                    max_cells  = 500L,
                    max_events = 200L,
                    na_colour  = "grey90",
                    title      = NULL, ...)
             standardGeneric("PlotHeatmap"))

#' Sashimi-style coverage plot for a splice event
#'
#' Draws junction arcs scaled by read count over a schematic gene structure.
#' Arcs are coloured by role: inclusion (blue) vs exclusion (red). Works in
#' both junction mode and transcript mode. Optionally faceted by a cell
#' metadata column.
#'
#' @param object A \code{MatisseObject}.
#' @param event_id Character. A single event ID; run
#'   \code{rownames(GetSeurat(obj)[["psi"]])} to list available IDs.
#' @param ... Additional arguments (see \code{\link{PlotSashimi}}).
#' @return A \code{ggplot} object.
#' @export
setGeneric("PlotSashimi",
           function(object, event_id,
                    cells     = NULL,
                    group_by  = NULL,
                    arc_scale = c("sqrt", "linear", "log"),
                    colours   = c(inclusion = "#4393c3", exclusion = "#d6604d"),
                    title     = NULL, ...)
             standardGeneric("PlotSashimi"))
