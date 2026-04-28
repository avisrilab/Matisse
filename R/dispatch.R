#' @include MatisseObject-class.R
#' @include MatisseObject-methods.R
NULL

# ---------------------------------------------------------------------------
# Hybrid method dispatch: Seurat and Signac generics on MatisseObject
#
# Each function is an S3 method that runs the equivalent Seurat/Signac
# function on the *embedded* Seurat object and wraps the result back into
# the MatisseObject transparently.
# ---------------------------------------------------------------------------

.seurat_forward <- function(FUN, object, ...) {
  result <- FUN(object@seurat, ...)
  if (inherits(result, "Seurat")) {
    object@seurat <- result
    return(object)
  }
  result
}

# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------

#' Normalise gene-expression counts for a MatisseObject
#'
#' Runs \code{\link[Seurat]{NormalizeData}} on the embedded Seurat object and
#' returns the updated \code{MatisseObject}. All splicing assays (junction,
#' psi, transcript) are unaffected.
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments forwarded to \code{\link[Seurat]{NormalizeData}}.
#' @return The updated \code{MatisseObject}.
#'
#' @seealso \code{\link{SCTransform.MatisseObject}}, \code{\link{ScaleData.MatisseObject}}
#' @rdname NormalizeData.MatisseObject
#' @importFrom Seurat NormalizeData
#' @method NormalizeData MatisseObject
#' @export
NormalizeData.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::NormalizeData, object, ...)
}

#' Scale gene-expression data for a MatisseObject
#'
#' Runs \code{\link[Seurat]{ScaleData}} on the embedded Seurat object and
#' returns the updated \code{MatisseObject}. Typically called after
#' \code{\link{NormalizeData.MatisseObject}} or
#' \code{\link{FindVariableFeatures.MatisseObject}}.
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments forwarded to \code{\link[Seurat]{ScaleData}}.
#' @return The updated \code{MatisseObject}.
#'
#' @seealso \code{\link{NormalizeData.MatisseObject}}, \code{\link{SCTransform.MatisseObject}}
#' @rdname ScaleData.MatisseObject
#' @importFrom Seurat ScaleData
#' @method ScaleData MatisseObject
#' @export
ScaleData.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::ScaleData, object, ...)
}

#' Identify highly variable features for a MatisseObject
#'
#' Runs \code{\link[Seurat]{FindVariableFeatures}} on the embedded Seurat
#' object and returns the updated \code{MatisseObject}. Identifies genes
#' whose expression varies most across cells — used to select features for
#' PCA.
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments forwarded to
#'   \code{\link[Seurat]{FindVariableFeatures}}.
#' @return The updated \code{MatisseObject}.
#'
#' @seealso \code{\link{RunPCA.MatisseObject}}, \code{\link{NormalizeData.MatisseObject}}
#' @rdname FindVariableFeatures.MatisseObject
#' @importFrom Seurat FindVariableFeatures
#' @method FindVariableFeatures MatisseObject
#' @export
FindVariableFeatures.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::FindVariableFeatures, object, ...)
}

#' SCTransform normalisation for MatisseObjects
#'
#' Runs \code{\link[Seurat]{SCTransform}} with mode-aware defaults.
#' In \strong{event mode} (long-read), normalises the \code{"transcript"}
#' assay by default, so that transcript-level abundances are variance-stabilised
#' before dimensionality reduction. In \strong{junction mode} (short-read),
#' normalises the active default assay (typically \code{"RNA"}).
#'
#' Override the target assay with the \code{assay} argument. PCA is
#' \strong{not} run automatically — call \code{\link{RunPCA.MatisseObject}}
#' after normalisation:
#'
#' \preformatted{
#' obj <- SCTransform(obj)                        # normalise
#' obj <- RunPCA(obj, assay = "SCT", npcs = 50)  # reduce
#' obj <- RunUMAP(obj, dims = 1:50)              # embed
#' }
#'
#' @param object A \code{MatisseObject}.
#' @param assay Character. Assay to normalise. Default: \code{"transcript"}
#'   in event mode; the active default assay in junction mode.
#' @param vars_to_regress Character vector. Covariates to regress out
#'   (e.g. \code{"percent.mt"}). Default: \code{NULL}.
#' @param verbose Logical. Print progress messages. Default: \code{TRUE}.
#' @param ... Additional arguments forwarded to \code{\link[Seurat]{SCTransform}}.
#' @return The updated \code{MatisseObject} with a new \code{"SCT"} assay.
#'
#' @seealso \code{\link{RunPCA.MatisseObject}}, \code{\link{NormalizeData.MatisseObject}}
#' @importFrom Seurat SCTransform
#' @method SCTransform MatisseObject
#' @export
SCTransform.MatisseObject <- function(object,
                                       assay           = NULL,
                                       vars_to_regress = NULL,
                                       verbose         = TRUE,
                                       ...) {
  if (is.null(assay)) {
    assay <- if (object@mode == "event") "transcript"
             else SeuratObject::DefaultAssay(object@seurat)
  }

  seu <- object@seurat
  if (is.null(.get_assay_safe(seu, assay))) {
    rlang::abort(paste0(
      "No '", assay, "' assay found. ",
      if (object@mode == "event")
        "Run CreateMatisseObject(transcript_counts=..., ioe_files=...) first."
      else
        "Run CreateMatisseObject(junction_counts=...) first, or pass `assay` explicitly."
    ))
  }

  SeuratObject::DefaultAssay(seu) <- assay

  if (verbose) {
    cli::cli_alert_info(
      "Running SCTransform on '{assay}' assay ({nrow(seu[[assay]])} features)...")
  }

  seu <- Seurat::SCTransform(
    seu,
    assay           = assay,
    vars.to.regress = vars_to_regress,
    verbose         = verbose,
    ...
  )

  object@seurat <- seu

  if (verbose) {
    cli::cli_alert_success("SCTransform complete. Run RunPCA() next.")
  }

  object
}

# ---------------------------------------------------------------------------
# Dimensionality reduction
# ---------------------------------------------------------------------------

#' Run PCA on a MatisseObject
#'
#' Runs \code{\link[Seurat]{RunPCA}} on the embedded Seurat object and
#' returns the updated \code{MatisseObject}. The PCA result is stored inside
#' the Seurat object and accessible via \code{GetSeurat(obj)}.
#'
#' Typical usage after \code{\link{SCTransform.MatisseObject}}:
#' \preformatted{
#' obj <- RunPCA(obj, assay = "SCT", npcs = 50)
#' }
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments forwarded to \code{\link[Seurat]{RunPCA}}
#'   (e.g. \code{assay}, \code{npcs}, \code{features}).
#' @return The updated \code{MatisseObject} with a \code{"pca"} reduction.
#'
#' @seealso \code{\link{RunUMAP.MatisseObject}}, \code{\link{SCTransform.MatisseObject}}
#' @rdname RunPCA.MatisseObject
#' @importFrom Seurat RunPCA
#' @method RunPCA MatisseObject
#' @export
RunPCA.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::RunPCA, object, ...)
}

#' Run UMAP on a MatisseObject
#'
#' Runs \code{\link[Seurat]{RunUMAP}} on the embedded Seurat object and
#' returns the updated \code{MatisseObject}. Call after
#' \code{\link{RunPCA.MatisseObject}} (or \code{\link{RunSVD.MatisseObject}}
#' for ATAC data). The resulting embedding is accessible via
#' \code{GetSeurat(obj)} and used by \code{\link{PlotUMAP}}.
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments forwarded to \code{\link[Seurat]{RunUMAP}}
#'   (e.g. \code{dims}, \code{reduction}).
#' @return The updated \code{MatisseObject} with a \code{"umap"} reduction.
#'
#' @seealso \code{\link{RunPCA.MatisseObject}}, \code{\link{PlotUMAP}}
#' @rdname RunUMAP.MatisseObject
#' @importFrom Seurat RunUMAP
#' @method RunUMAP MatisseObject
#' @export
RunUMAP.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::RunUMAP, object, ...)
}

#' Run t-SNE on a MatisseObject
#'
#' Runs \code{\link[Seurat]{RunTSNE}} on the embedded Seurat object and
#' returns the updated \code{MatisseObject}. An alternative to
#' \code{\link{RunUMAP.MatisseObject}} for 2-D cell embedding.
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments forwarded to \code{\link[Seurat]{RunTSNE}}
#'   (e.g. \code{dims}, \code{perplexity}).
#' @return The updated \code{MatisseObject} with a \code{"tsne"} reduction.
#'
#' @seealso \code{\link{RunUMAP.MatisseObject}}, \code{\link{RunPCA.MatisseObject}}
#' @rdname RunTSNE.MatisseObject
#' @importFrom Seurat RunTSNE
#' @method RunTSNE MatisseObject
#' @export
RunTSNE.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::RunTSNE, object, ...)
}

# ---------------------------------------------------------------------------
# Clustering
# ---------------------------------------------------------------------------

#' Compute a shared nearest-neighbour graph for a MatisseObject
#'
#' Runs \code{\link[Seurat]{FindNeighbors}} on the embedded Seurat object and
#' returns the updated \code{MatisseObject}. Typically called after
#' \code{\link{RunPCA.MatisseObject}} and before
#' \code{\link{FindClusters.MatisseObject}}.
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments forwarded to
#'   \code{\link[Seurat]{FindNeighbors}} (e.g. \code{dims}, \code{reduction}).
#' @return The updated \code{MatisseObject} with neighbour graph stored
#'   inside the embedded Seurat object.
#'
#' @seealso \code{\link{FindClusters.MatisseObject}}, \code{\link{RunPCA.MatisseObject}}
#' @rdname FindNeighbors.MatisseObject
#' @importFrom Seurat FindNeighbors
#' @method FindNeighbors MatisseObject
#' @export
FindNeighbors.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::FindNeighbors, object, ...)
}

#' Cluster cells in a MatisseObject
#'
#' Runs \code{\link[Seurat]{FindClusters}} on the embedded Seurat object and
#' returns the updated \code{MatisseObject}. Cluster assignments are stored
#' in \code{MatisseMeta(obj)$seurat_clusters} and are immediately available
#' for \code{\link{PlotUMAP}} and \code{\link{PlotViolin}}.
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments forwarded to
#'   \code{\link[Seurat]{FindClusters}} (e.g. \code{resolution}).
#' @return The updated \code{MatisseObject} with \code{seurat_clusters}
#'   added to cell metadata.
#'
#' @seealso \code{\link{FindNeighbors.MatisseObject}}, \code{\link{PlotUMAP}}
#' @rdname FindClusters.MatisseObject
#' @importFrom Seurat FindClusters
#' @method FindClusters MatisseObject
#' @export
FindClusters.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::FindClusters, object, ...)
}

#' Find differentially expressed markers for a MatisseObject
#'
#' Runs \code{\link[Seurat]{FindMarkers}} on the embedded Seurat object.
#' Unlike most dispatch methods, this returns a \code{data.frame} of
#' marker statistics rather than an updated \code{MatisseObject}.
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments forwarded to
#'   \code{\link[Seurat]{FindMarkers}} (e.g. \code{ident.1}, \code{ident.2},
#'   \code{group.by}, \code{features}).
#' @return A \code{data.frame} of marker genes with columns
#'   \code{p_val}, \code{avg_log2FC}, \code{pct.1}, \code{pct.2},
#'   \code{p_val_adj}.
#'
#' @seealso \code{\link{FindClusters.MatisseObject}}
#' @rdname FindMarkers.MatisseObject
#' @importFrom Seurat FindMarkers
#' @method FindMarkers MatisseObject
#' @export
FindMarkers.MatisseObject <- function(object, ...) {
  Seurat::FindMarkers(object@seurat, ...)
}

#' Add metadata columns to a MatisseObject
#'
#' Runs \code{\link[SeuratObject]{AddMetaData}} on the embedded Seurat object
#' and returns the updated \code{MatisseObject}. New columns are immediately
#' accessible via \code{\link{MatisseMeta}} and the \code{$} operator.
#'
#' @param object A \code{MatisseObject}.
#' @param metadata A named vector, \code{data.frame}, or list of columns to add.
#' @param col.name Character. Column name to use when \code{metadata} is a
#'   vector. Default: \code{NULL}.
#' @param ... Additional arguments forwarded to
#'   \code{\link[SeuratObject]{AddMetaData}}.
#' @return The updated \code{MatisseObject} with new metadata columns.
#'
#' @seealso \code{\link{MatisseMeta}}, \code{\link{AddIsoformMetadata}}
#' @rdname AddMetaData.MatisseObject
#' @importFrom SeuratObject AddMetaData
#' @method AddMetaData MatisseObject
#' @export
AddMetaData.MatisseObject <- function(object, metadata, col.name = NULL, ...) {
  result <- SeuratObject::AddMetaData(object@seurat, metadata = metadata,
                                       col.name = col.name, ...)
  if (inherits(result, "Seurat")) {
    object@seurat <- result
    return(object)
  }
  result
}

# ---------------------------------------------------------------------------
# Signac methods (ATAC / multiome)
# ---------------------------------------------------------------------------

#' Run TF-IDF normalisation for a MatisseObject
#'
#' Runs \code{\link[Signac]{RunTFIDF}} on the embedded Seurat object and
#' returns the updated \code{MatisseObject}. Used for ATAC-seq peak counts
#' in multiome datasets before \code{\link{RunSVD.MatisseObject}}.
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments forwarded to
#'   \code{\link[Signac]{RunTFIDF}} (e.g. \code{assay}, \code{method}).
#' @return The updated \code{MatisseObject} with TF-IDF normalised counts.
#'
#' @seealso \code{\link{RunSVD.MatisseObject}}, \code{\link{FindTopFeatures.MatisseObject}}
#' @rdname RunTFIDF.MatisseObject
#' @importFrom Signac RunTFIDF
#' @method RunTFIDF MatisseObject
#' @export
RunTFIDF.MatisseObject <- function(object, ...) {
  .seurat_forward(Signac::RunTFIDF, object, ...)
}

#' Run SVD (LSI) for a MatisseObject
#'
#' Runs \code{\link[Signac]{RunSVD}} on the embedded Seurat object and
#' returns the updated \code{MatisseObject}. Performs Latent Semantic
#' Indexing (LSI) — the standard dimensionality-reduction step for ATAC-seq
#' data. Typically called after \code{\link{RunTFIDF.MatisseObject}}.
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments forwarded to \code{\link[Signac]{RunSVD}}
#'   (e.g. \code{n}, \code{assay}).
#' @return The updated \code{MatisseObject} with an \code{"lsi"} reduction.
#'
#' @seealso \code{\link{RunTFIDF.MatisseObject}}, \code{\link{RunUMAP.MatisseObject}}
#' @rdname RunSVD.MatisseObject
#' @importFrom Signac RunSVD
#' @method RunSVD MatisseObject
#' @export
RunSVD.MatisseObject <- function(object, ...) {
  .seurat_forward(Signac::RunSVD, object, ...)
}

#' Find highly variable ATAC-seq features for a MatisseObject
#'
#' Runs \code{\link[Signac]{FindTopFeatures}} on the embedded Seurat object
#' and returns the updated \code{MatisseObject}. Selects the most accessible
#' peaks for downstream LSI / UMAP.
#'
#' @param object A \code{MatisseObject}.
#' @param ... Additional arguments forwarded to
#'   \code{\link[Signac]{FindTopFeatures}} (e.g. \code{min.cutoff}).
#' @return The updated \code{MatisseObject} with top ATAC features flagged.
#'
#' @seealso \code{\link{RunTFIDF.MatisseObject}}, \code{\link{RunSVD.MatisseObject}}
#' @rdname FindTopFeatures.MatisseObject
#' @importFrom Signac FindTopFeatures
#' @method FindTopFeatures MatisseObject
#' @export
FindTopFeatures.MatisseObject <- function(object, ...) {
  .seurat_forward(Signac::FindTopFeatures, object, ...)
}
