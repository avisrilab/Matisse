#' @include MatisseObject-class.R
#' @include MatisseObject-methods.R
NULL

# ---------------------------------------------------------------------------
# Hybrid method dispatch: Seurat and Signac generics on MatisseObject
#
# When a Seurat or Signac function is called directly on a MatisseObject
# (e.g. RunPCA(matisse_obj)), these S3 forwarding methods intercept the call,
# run the function on the embedded Seurat object, and—when the result is a
# Seurat object—wrap it back into a MatisseObject transparently.
#
# Functions that return something other than a Seurat object (e.g. Embeddings,
# FindMarkers) return their result directly.
#
# The $ operator on MatisseObject provides the same forwarding for any
# Seurat/Signac function not listed here (see MatisseObject-methods.R).
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
# SeuratObject generics (always available; SeuratObject is in Imports)
# ---------------------------------------------------------------------------

#' @importFrom SeuratObject NormalizeData
#' @export
NormalizeData.MatisseObject <- function(object, ...) {
  .seurat_forward(SeuratObject::NormalizeData, object, ...)
}

#' @importFrom SeuratObject ScaleData
#' @export
ScaleData.MatisseObject <- function(object, ...) {
  .seurat_forward(SeuratObject::ScaleData, object, ...)
}

#' @importFrom SeuratObject FindVariableFeatures
#' @export
FindVariableFeatures.MatisseObject <- function(object, ...) {
  .seurat_forward(SeuratObject::FindVariableFeatures, object, ...)
}

#' @importFrom SeuratObject AddMetaData
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
# Seurat generics (Seurat is in Imports)
# ---------------------------------------------------------------------------

#' @importFrom Seurat RunPCA
#' @export
RunPCA.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::RunPCA, object, ...)
}

#' @importFrom Seurat RunUMAP
#' @export
RunUMAP.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::RunUMAP, object, ...)
}

#' @importFrom Seurat RunTSNE
#' @export
RunTSNE.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::RunTSNE, object, ...)
}

#' @importFrom Seurat FindNeighbors
#' @export
FindNeighbors.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::FindNeighbors, object, ...)
}

#' @importFrom Seurat FindClusters
#' @export
FindClusters.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::FindClusters, object, ...)
}

#' @importFrom Seurat SCTransform
#' @export
SCTransform.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::SCTransform, object, ...)
}

#' @importFrom Seurat FindMarkers
#' @export
FindMarkers.MatisseObject <- function(object, ...) {
  Seurat::FindMarkers(object@seurat, ...)
}

# ---------------------------------------------------------------------------
# Signac generics (Signac is in Imports)
# ---------------------------------------------------------------------------

#' @importFrom Signac RunTFIDF
#' @export
RunTFIDF.MatisseObject <- function(object, ...) {
  .seurat_forward(Signac::RunTFIDF, object, ...)
}

#' @importFrom Signac RunSVD
#' @export
RunSVD.MatisseObject <- function(object, ...) {
  .seurat_forward(Signac::RunSVD, object, ...)
}

#' @importFrom Signac FindTopFeatures
#' @export
FindTopFeatures.MatisseObject <- function(object, ...) {
  .seurat_forward(Signac::FindTopFeatures, object, ...)
}
