#' @include MatisseObject-class.R
#' @include MatisseObject-methods.R
NULL

# ---------------------------------------------------------------------------
# Hybrid method dispatch: Seurat and Signac generics on MatisseObject
#
# Seurat/Signac generics are registered as both S3 methods (via @method) and
# S4 methods (via setMethod) so that direct calls like RunPCA(obj, ...) work
# regardless of which dispatch path is taken.  When the result is a Seurat
# object it is wrapped back into the MatisseObject transparently.
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

#' @importFrom SeuratObject NormalizeData
#' @method NormalizeData MatisseObject
#' @export
NormalizeData.MatisseObject <- function(object, ...) {
  .seurat_forward(SeuratObject::NormalizeData, object, ...)
}
setMethod("NormalizeData", "MatisseObject", NormalizeData.MatisseObject)

#' @importFrom SeuratObject ScaleData
#' @method ScaleData MatisseObject
#' @export
ScaleData.MatisseObject <- function(object, ...) {
  .seurat_forward(SeuratObject::ScaleData, object, ...)
}
setMethod("ScaleData", "MatisseObject", ScaleData.MatisseObject)

#' @importFrom SeuratObject FindVariableFeatures
#' @method FindVariableFeatures MatisseObject
#' @export
FindVariableFeatures.MatisseObject <- function(object, ...) {
  .seurat_forward(SeuratObject::FindVariableFeatures, object, ...)
}
setMethod("FindVariableFeatures", "MatisseObject", FindVariableFeatures.MatisseObject)

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
setMethod("AddMetaData", "MatisseObject", AddMetaData.MatisseObject)

# ---------------------------------------------------------------------------
# Seurat generics
# ---------------------------------------------------------------------------

#' @importFrom Seurat RunPCA
#' @method RunPCA MatisseObject
#' @export
RunPCA.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::RunPCA, object, ...)
}
setMethod("RunPCA", "MatisseObject", RunPCA.MatisseObject)

#' @importFrom Seurat RunUMAP
#' @method RunUMAP MatisseObject
#' @export
RunUMAP.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::RunUMAP, object, ...)
}
setMethod("RunUMAP", "MatisseObject", RunUMAP.MatisseObject)

#' @importFrom Seurat RunTSNE
#' @method RunTSNE MatisseObject
#' @export
RunTSNE.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::RunTSNE, object, ...)
}
setMethod("RunTSNE", "MatisseObject", RunTSNE.MatisseObject)

#' @importFrom Seurat FindNeighbors
#' @method FindNeighbors MatisseObject
#' @export
FindNeighbors.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::FindNeighbors, object, ...)
}
setMethod("FindNeighbors", "MatisseObject", FindNeighbors.MatisseObject)

#' @importFrom Seurat FindClusters
#' @method FindClusters MatisseObject
#' @export
FindClusters.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::FindClusters, object, ...)
}
setMethod("FindClusters", "MatisseObject", FindClusters.MatisseObject)

#' @importFrom Seurat SCTransform
#' @method SCTransform MatisseObject
#' @export
SCTransform.MatisseObject <- function(object, ...) {
  .seurat_forward(Seurat::SCTransform, object, ...)
}
setMethod("SCTransform", "MatisseObject", SCTransform.MatisseObject)

#' @importFrom Seurat FindMarkers
#' @method FindMarkers MatisseObject
#' @export
FindMarkers.MatisseObject <- function(object, ...) {
  Seurat::FindMarkers(object@seurat, ...)
}
setMethod("FindMarkers", "MatisseObject", FindMarkers.MatisseObject)

# ---------------------------------------------------------------------------
# Signac generics
# ---------------------------------------------------------------------------

#' @importFrom Signac RunTFIDF
#' @method RunTFIDF MatisseObject
#' @export
RunTFIDF.MatisseObject <- function(object, ...) {
  .seurat_forward(Signac::RunTFIDF, object, ...)
}
setMethod("RunTFIDF", "MatisseObject", RunTFIDF.MatisseObject)

#' @importFrom Signac RunSVD
#' @method RunSVD MatisseObject
#' @export
RunSVD.MatisseObject <- function(object, ...) {
  .seurat_forward(Signac::RunSVD, object, ...)
}
setMethod("RunSVD", "MatisseObject", RunSVD.MatisseObject)

#' @importFrom Signac FindTopFeatures
#' @method FindTopFeatures MatisseObject
#' @export
FindTopFeatures.MatisseObject <- function(object, ...) {
  .seurat_forward(Signac::FindTopFeatures, object, ...)
}
setMethod("FindTopFeatures", "MatisseObject", FindTopFeatures.MatisseObject)
