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
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments (unused).
#' @return A sparse matrix (cells x events) of PSI values.
#' @export
setGeneric("GetPSI", function(object, ...) standardGeneric("GetPSI"))

#' Set the PSI matrix
#' @param object A \code{MatisseObject}.
#' @param value A sparse matrix (cells x events) of PSI values.
#' @return The updated \code{MatisseObject}.
#' @export
setGeneric("SetPSI", function(object, value) standardGeneric("SetPSI"))

#' Get raw junction count matrix
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments (unused).
#' @return A sparse matrix (cells x junctions) of read counts.
#' @export
setGeneric("GetJunctionCounts",
           function(object, ...) standardGeneric("GetJunctionCounts"))

#' Get inclusion read count matrix
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments (unused).
#' @return A sparse matrix (cells x events) of inclusion read counts.
#' @export
setGeneric("GetInclusionCounts",
           function(object, ...) standardGeneric("GetInclusionCounts"))

#' Get exclusion read count matrix
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments (unused).
#' @return A sparse matrix (cells x events) of exclusion read counts.
#' @export
setGeneric("GetExclusionCounts",
           function(object, ...) standardGeneric("GetExclusionCounts"))

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

#' Get or set cell-level isoform metadata
#'
#' @param object A \code{MatisseObject}.
#' @param value A \code{data.frame} of cell-level metadata to assign (for the
#'   setter). Rownames must match cell barcodes.
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

#' Add columns to the isoform metadata
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
#' @param object A \code{MatisseObject} or a sparse count matrix
#'   (cells x junctions).
#' @param events A \code{data.frame} defining splice events (see
#'   \code{\link{CalculatePSI}} for required columns).
#' @param ... Additional arguments passed to the method.
#' @return A \code{MatisseObject} (when given one) or a PSI matrix.
#' @export
setGeneric("CalculatePSI",
           function(object, events = NULL, ...) standardGeneric("CalculatePSI"))

#' Compute per-cell isoform QC metrics
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments (unused).
#' @return The updated \code{MatisseObject} with QC columns added to
#'   isoform metadata.
#' @export
setGeneric("ComputeIsoformQC",
           function(object, ...) standardGeneric("ComputeIsoformQC"))

#' Filter cells by isoform QC thresholds
#'
#' @param object A \code{MatisseObject}.
#' @param ... Named numeric thresholds; see \code{\link{FilterCells}}.
#' @return The filtered \code{MatisseObject}.
#' @export
setGeneric("FilterCells",
           function(object, ...) standardGeneric("FilterCells"))

#' Filter splice events by coverage or variance
#'
#' @param object A \code{MatisseObject}.
#' @param ... Named thresholds; see \code{\link{FilterEvents}}.
#' @return The filtered \code{MatisseObject}.
#' @export
setGeneric("FilterEvents",
           function(object, ...) standardGeneric("FilterEvents"))

# ---------------------------------------------------------------------------
# Visualization generics
# ---------------------------------------------------------------------------

#' Heatmap of PSI values across cells and events
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments (see \code{\link{PlotPSIHeatmap}}).
#' @return A \code{ggplot} or \code{pheatmap} object.
#' @export
setGeneric("PlotPSIHeatmap",
           function(object, ...) standardGeneric("PlotPSIHeatmap"))

#' UMAP plot colored by PSI of a specific splice event
#' @param object A \code{MatisseObject}.
#' @param event_id Character. Name of the splice event to visualize.
#' @param ... Additional arguments (see \code{\link{PlotPSIUMAP}}).
#' @return A \code{ggplot} object.
#' @export
setGeneric("PlotPSIUMAP",
           function(object, event_id, ...) standardGeneric("PlotPSIUMAP"))

#' Violin plot of PSI values split by cell group
#' @param object A \code{MatisseObject}.
#' @param event_id Character. Name of the splice event to visualize.
#' @param ... Additional arguments (see \code{\link{PlotPSIViolin}}).
#' @return A \code{ggplot} object.
#' @export
setGeneric("PlotPSIViolin",
           function(object, event_id, ...) standardGeneric("PlotPSIViolin"))

#' Per-cell junction coverage plot for a gene
#' @param object A \code{MatisseObject}.
#' @param gene Character. Gene name.
#' @param ... Additional arguments (see \code{\link{PlotJunctionCoverage}}).
#' @return A \code{ggplot} object.
#' @export
setGeneric("PlotJunctionCoverage",
           function(object, gene, ...) standardGeneric("PlotJunctionCoverage"))

#' Violin/ridge plot of isoform QC metrics
#' @param object A \code{MatisseObject}.
#' @param features Character vector of QC metric names to plot.
#' @param ... Additional arguments (see \code{\link{PlotQCMetrics}}).
#' @return A \code{ggplot} object.
#' @export
setGeneric("PlotQCMetrics",
           function(object, features = NULL, ...) standardGeneric("PlotQCMetrics"))
