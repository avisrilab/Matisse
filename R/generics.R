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
#' @return A sparse matrix (cells × events) of PSI values in \eqn{[0,1]}.
#'   \code{NULL} if no \code{"psi"} assay exists yet.
#' @export
setGeneric("GetPSI", function(object, ...) standardGeneric("GetPSI"))

#' Set the PSI matrix
#'
#' Replaces the \code{"data"} layer of the \code{"psi"} \code{Assay5}
#' inside the embedded Seurat object.
#'
#' @param object A \code{MatisseObject}.
#' @param value A sparse matrix (cells × events) of PSI values.
#' @return The updated \code{MatisseObject}.
#' @export
setGeneric("SetPSI", function(object, value) standardGeneric("SetPSI"))

#' Get raw junction count matrix
#'
#' Retrieves the per-junction read counts from the \code{"junction"}
#' \code{Assay5} stored inside the embedded Seurat object.
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments (unused).
#' @return A sparse matrix (cells × junctions) of read counts, or \code{NULL}
#'   if the object is in event mode or no junction assay exists.
#' @export
setGeneric("GetJunctionCounts",
           function(object, ...) standardGeneric("GetJunctionCounts"))

#' Get inclusion read count matrix
#'
#' Retrieves inclusion counts from the \code{"counts"} layer of the
#' \code{"psi"} \code{Assay5}.
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments (unused).
#' @return A sparse matrix (cells × events) of inclusion read counts.
#' @export
setGeneric("GetInclusionCounts",
           function(object, ...) standardGeneric("GetInclusionCounts"))

#' Get exclusion read count matrix
#'
#' Retrieves exclusion counts from the \code{"exclusion"} layer of the
#' \code{"psi"} \code{Assay5}.
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments (unused).
#' @return A sparse matrix (cells × events) of exclusion read counts.
#' @export
setGeneric("GetExclusionCounts",
           function(object, ...) standardGeneric("GetExclusionCounts"))

#' Get transcript count matrix
#'
#' Retrieves raw transcript counts from the \code{"transcript"} \code{Assay5}
#' stored inside the embedded Seurat object.
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments (unused).
#' @return A sparse matrix (transcripts × cells) of raw counts, or \code{NULL}
#'   if no \code{"transcript"} assay exists.
#' @export
setGeneric("GetTranscriptCounts",
           function(object, ...) standardGeneric("GetTranscriptCounts"))

#' Get splice event annotation table
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments (unused).
#' @return A \code{data.frame} of splice event metadata.
#' @export
setGeneric("GetEventData",
           function(object, ...) standardGeneric("GetEventData"))

#' Get junction annotation table
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments (unused).
#' @return A \code{data.frame} of junction metadata.
#' @export
setGeneric("GetJunctionData",
           function(object, ...) standardGeneric("GetJunctionData"))

#' Get or set cell-level metadata
#'
#' Returns the full \code{meta.data} of the embedded Seurat object, which
#' includes all per-cell QC metrics and annotations added by Matisse (e.g.
#' \code{n_junctions_detected}, \code{mean_psi}) alongside standard Seurat
#' columns. Use \code{AddIsoformMetadata()} to add new columns.
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

#' Add columns to the cell metadata
#'
#' Adds new columns to the embedded Seurat object's \code{meta.data}. This is
#' the standard way to attach per-cell isoform QC or other annotations to a
#' \code{MatisseObject}.
#'
#' @param object A \code{MatisseObject}.
#' @param metadata A named \code{data.frame} or named numeric/character vector.
#'   Rownames (or names) must match cell barcodes.
#' @param ... Additional arguments (unused).
#' @return The updated \code{MatisseObject}.
#' @export
setGeneric("AddIsoformMetadata",
           function(object, metadata, ...) standardGeneric("AddIsoformMetadata"))

# ---------------------------------------------------------------------------
# Analysis generics
# ---------------------------------------------------------------------------

#' Calculate PSI matrix from junction counts
#'
#' @param object A \code{MatisseObject} in junction mode, or a sparse count
#'   matrix (cells × junctions).
#' @param events A \code{data.frame} defining splice events (see
#'   \code{\link{CalculatePSI}} for required columns).
#' @param ... Additional arguments passed to the method.
#' @return A \code{MatisseObject} (when given one) or a PSI matrix.
#' @export
setGeneric("CalculatePSI",
           function(object, events = NULL, min_coverage = 5L,
                    na_fill = NA_real_, verbose = TRUE, ...)
             standardGeneric("CalculatePSI"))

#' Compute per-cell isoform QC metrics
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments (unused).
#' @return The updated \code{MatisseObject} with QC columns added to
#'   cell metadata.
#' @export
setGeneric("ComputeIsoformQC",
           function(object, min_coverage = 5L, verbose = TRUE, ...)
             standardGeneric("ComputeIsoformQC"))

#' Filter cells by isoform QC thresholds
#'
#' @param object A \code{MatisseObject}.
#' @param min_junctions Integer. Minimum \code{n_junctions_detected}. Default: \code{NULL}.
#' @param max_junctions Integer. Maximum \code{n_junctions_detected}. Default: \code{NULL}.
#' @param min_junction_reads Integer. Minimum \code{total_junction_reads}. Default: \code{NULL}.
#' @param max_junction_reads Integer. Maximum \code{total_junction_reads}. Default: \code{NULL}.
#' @param min_pct_covered Numeric (0–100). Minimum \code{pct_events_covered}. Default: \code{NULL}.
#' @param custom_filters Named list of \code{c(min, max)} bounds for arbitrary metadata columns.
#' @param verbose Logical. Default: \code{TRUE}.
#' @return The filtered \code{MatisseObject}.
#' @export
setGeneric("FilterCells",
           function(object,
                    min_junctions      = NULL,
                    max_junctions      = NULL,
                    min_junction_reads = NULL,
                    max_junction_reads = NULL,
                    min_pct_covered    = NULL,
                    custom_filters     = NULL,
                    verbose            = TRUE, ...)
             standardGeneric("FilterCells"))

#' Filter splice events by coverage or variance
#'
#' @param object A \code{MatisseObject}.
#' @param min_cells_covered Integer. Minimum cells with non-NA PSI. Default: \code{10}.
#' @param min_psi_variance Numeric. Minimum PSI variance across covered cells. Default: \code{NULL}.
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
#' Overlays the value of a feature (PSI event, junction count, or gene
#' expression) on the UMAP embedding stored in the embedded Seurat object.
#' Pass an event ID for PSI, a junction ID for junction counts, or a gene
#' name for expression.
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
#' @param object A \code{MatisseObject}.
#' @param feature Character. Feature to visualise (PSI event, junction, or gene).
#' @param group_by Character. Metadata column to split cells by. Default: \code{"seurat_clusters"}.
#' @param colours Named character vector mapping group levels to colours. Default: \code{NULL}.
#' @param add_points Logical. Overlay jittered cell values. Default: \code{FALSE}.
#' @param title Character. Plot title. Default: feature name.
#' @return A \code{ggplot} object.
#' @export
setGeneric("PlotViolin",
           function(object, feature,
                    group_by   = "seurat_clusters",
                    colours    = NULL,
                    add_points = FALSE,
                    title      = NULL, ...)
             standardGeneric("PlotViolin"))

#' Heatmap of feature values across cells and events
#'
#' @param object A \code{MatisseObject}.
#' @param events Character vector of event IDs. Default: all events.
#' @param cells Character vector of cell barcodes. Default: all cells.
#' @param group_by Character. Metadata column to annotate and order cells. Default: \code{NULL}.
#' @param max_cells Integer. Downsample cap before plotting. Default: \code{500}.
#' @param na_colour Character. Colour for \code{NA} entries. Default: \code{"grey90"}.
#' @return A \code{ggplot} object.
#' @export
setGeneric("PlotHeatmap",
           function(object,
                    events    = NULL,
                    cells     = NULL,
                    group_by  = NULL,
                    max_cells = 500L,
                    na_colour = "grey90", ...)
             standardGeneric("PlotHeatmap"))

#' Junction coverage bar plot for a gene
#'
#' @param object A \code{MatisseObject} in junction mode.
#' @param gene Character. Gene name.
#' @param ... Additional arguments (see \code{\link{PlotCoverage}}).
#' @return A \code{ggplot} object.
#' @export
setGeneric("PlotCoverage",
           function(object, gene, cells = NULL, log_scale = FALSE, ...)
             standardGeneric("PlotCoverage"))

#' Violin/ridge plot of isoform QC metrics
#' @param object A \code{MatisseObject}.
#' @param features Character vector of QC metric names to plot.
#' @param group_by Character. Metadata column to split cells by. Default: \code{NULL}.
#' @param ncol Integer. Number of facet columns. Default: \code{2}.
#' @param ... Additional arguments passed to methods.
#' @return A \code{ggplot} object.
#' @export
setGeneric("PlotQCMetrics",
           function(object, features = NULL, group_by = NULL, ncol = 2L, ...)
             standardGeneric("PlotQCMetrics"))

#' Sashimi-style coverage plot for a splice event
#'
#' Draws junction arcs scaled by read count over a schematic gene structure.
#' Arcs are coloured by role: inclusion (blue) vs exclusion (red). Works in
#' both junction mode (per-junction counts) and event mode (aggregated
#' inclusion/exclusion counts). Optionally faceted by a cell metadata column.
#'
#' @param object A \code{MatisseObject}.
#' @param event_id Character. Event ID as stored in \code{event_data}.
#' @param ... Additional arguments (see \code{\link{CoveragePlot}}).
#' @return A \code{ggplot} object.
#' @export
setGeneric("CoveragePlot",
           function(object, event_id,
                    cells     = NULL,
                    group_by  = NULL,
                    arc_scale = c("sqrt", "linear", "log"),
                    colours   = c(inclusion = "#4393c3", exclusion = "#d6604d"),
                    title     = NULL, ...)
             standardGeneric("CoveragePlot"))
