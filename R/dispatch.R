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

#' @importFrom Seurat SCTransform
#' @method SCTransform MatisseObject
#' @export
SCTransform.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::SCTransform, object, ...)
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
