#' @include MatisseObject-class.R
#' @include MatisseObject-methods.R
NULL

# ---------------------------------------------------------------------------
# Hybrid method dispatch: Seurat and Signac generics on MatisseObject
#
# Each function is registered as an S3 method via @method so that direct
# calls like RunPCA(obj, ...) dispatch here automatically.  When the result
# is a Seurat object it is wrapped back into the MatisseObject transparently.
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
# SeuratObject generics
# ---------------------------------------------------------------------------

#' @importFrom Seurat NormalizeData
#' @method NormalizeData MatisseObject
#' @export
NormalizeData.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::NormalizeData, object, ...)
}

#' @importFrom Seurat ScaleData
#' @method ScaleData MatisseObject
#' @export
ScaleData.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::ScaleData, object, ...)
}

#' @importFrom Seurat FindVariableFeatures
#' @method FindVariableFeatures MatisseObject
#' @export
FindVariableFeatures.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::FindVariableFeatures, object, ...)
}

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
# Seurat generics
# ---------------------------------------------------------------------------

#' @importFrom Seurat RunPCA
#' @method RunPCA MatisseObject
#' @export
RunPCA.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::RunPCA, object, ...)
}

#' @importFrom Seurat RunUMAP
#' @method RunUMAP MatisseObject
#' @export
RunUMAP.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::RunUMAP, object, ...)
}

#' @importFrom Seurat RunTSNE
#' @method RunTSNE MatisseObject
#' @export
RunTSNE.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::RunTSNE, object, ...)
}

#' @importFrom Seurat FindNeighbors
#' @method FindNeighbors MatisseObject
#' @export
FindNeighbors.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::FindNeighbors, object, ...)
}

#' @importFrom Seurat FindClusters
#' @method FindClusters MatisseObject
#' @export
FindClusters.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::FindClusters, object, ...)
}

#' SCTransform normalisation for MatisseObjects
#'
#' Runs \code{\link[Seurat]{SCTransform}} with mode-aware defaults. In
#' \strong{event mode}, normalises the \code{"transcript"} assay. In
#' \strong{junction mode}, normalises the active default assay (usually
#' \code{"RNA"}). Override with the \code{assay} argument.
#'
#' PCA is \strong{not} run automatically. Call \code{RunPCA()} on the
#' resulting object after normalisation:
#' \preformatted{
#' obj <- SCTransform(obj)
#' obj <- RunPCA(obj, assay = "SCT", npcs = 50)
#' obj <- RunUMAP(obj, dims = 1:50)
#' }
#'
#' @param object A \code{MatisseObject}.
#' @param assay Character. Assay to normalise. Default: \code{"transcript"}
#'   in event mode; the active default assay in junction mode.
#' @param vars_to_regress Character vector. Variables to regress out.
#'   Default: \code{NULL}.
#' @param verbose Logical. Default: \code{TRUE}.
#' @param ... Additional arguments forwarded to \code{\link[Seurat]{SCTransform}}.
#' @return The updated \code{MatisseObject}.
#'
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

#' @importFrom Seurat FindMarkers
#' @method FindMarkers MatisseObject
#' @export
FindMarkers.MatisseObject <- function(object, ...) {
  Seurat::FindMarkers(object@seurat, ...)
}

# ---------------------------------------------------------------------------
# Signac generics
# ---------------------------------------------------------------------------

#' @importFrom Signac RunTFIDF
#' @method RunTFIDF MatisseObject
#' @export
RunTFIDF.MatisseObject <- function(object, ...) {
  .seurat_forward(Signac::RunTFIDF, object, ...)
}

#' @importFrom Signac RunSVD
#' @method RunSVD MatisseObject
#' @export
RunSVD.MatisseObject <- function(object, ...) {
  .seurat_forward(Signac::RunSVD, object, ...)
}

#' @importFrom Signac FindTopFeatures
#' @method FindTopFeatures MatisseObject
#' @export
FindTopFeatures.MatisseObject <- function(object, ...) {
  .seurat_forward(Signac::FindTopFeatures, object, ...)
}
