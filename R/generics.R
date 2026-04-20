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
#' \code{ChromatinAssay} stored inside the embedded Seurat object.
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments (unused).
#' @return A sparse matrix (cells × events) of PSI values in \eqn{[0,1]}.
#'   \code{NULL} if no \code{"psi"} assay exists yet.
#' @export
setGeneric("GetPSI", function(object, ...) standardGeneric("GetPSI"))

#' Set the PSI matrix
#'
#' Replaces the \code{"data"} layer of the \code{"psi"} \code{ChromatinAssay}
#' inside the embedded Seurat object.
#'
#' @param object A \code{MatisseObject}.
#' @param value A sparse matrix (cells × events) of PSI values.
#' @return The updated \code{MatisseObject}.
#' @export
setGeneric("SetPSI", function(object, value) standardGeneric("SetPSI"))

#' Get raw junction count matrix
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments (unused).
#' @return A sparse matrix (cells × junctions) of read counts.
#' @export
setGeneric("GetJunctionCounts",
           function(object, ...) standardGeneric("GetJunctionCounts"))

#' Get inclusion read count matrix
#'
#' Retrieves inclusion counts from the \code{"counts"} layer of the
#' \code{"psi"} \code{ChromatinAssay}.
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
#' \code{"psi"} \code{ChromatinAssay}.
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
#'   (cells × junctions).
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

#' SCTransform normalization for the transcript assay
#'
#' Runs \code{SCTransform} on the \code{"transcript"} assay of the embedded
#' Seurat object and follows it with \code{RunPCA} using a larger number of
#' principal components, giving better isoform-level cluster resolution.
#'
#' @param object A \code{MatisseObject} with a \code{"transcript"} assay.
#' @param ... Additional arguments passed to the method.
#' @return The updated \code{MatisseObject}.
#' @export
setGeneric("SCTransformTranscripts",
           function(object, ...) standardGeneric("SCTransformTranscripts"))

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
