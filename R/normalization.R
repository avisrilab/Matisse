#' @include MatisseObject-class.R
#' @include generics.R
NULL

# ---------------------------------------------------------------------------
# SCTransformTranscripts
# ---------------------------------------------------------------------------

#' SCTransform normalization for the transcript assay
#'
#' Runs \code{\link[Seurat]{SCTransform}} on the \code{"transcript"}
#' \code{Assay5} of the embedded Seurat object, followed by
#' \code{\link[Seurat]{RunPCA}} with \code{n_pca_dims} components.
#'
#' SCTransform is better suited than standard log-normalization for
#' transcript-level count data: it corrects for sequencing depth while
#' preserving biological variance. Using a higher number of PCA dimensions
#' (default 50) has been found to resolve more isoform-level clusters than the
#' standard 20.
#'
#' After this step you can call \code{\link[Seurat]{FindNeighbors}},
#' \code{\link[Seurat]{FindClusters}}, and \code{\link[Seurat]{RunUMAP}}
#' directly on the \code{MatisseObject} via the built-in method forwarding
#' (see \code{\link{MatisseObject-class}}).
#'
#' @param object A \code{\linkS4class{MatisseObject}} with a
#'   \code{"transcript"} assay. Use
#'   \code{\link{CreateMatisseObjectFromTranscripts}} or
#'   \code{CreateMatisseObject(transcript_counts = ...)} to create one.
#' @param n_pca_dims Integer. Number of principal components to compute.
#'   Default: \code{50}.
#' @param vars_to_regress Character vector of variables to regress out during
#'   SCTransform (e.g. \code{"percent.mt"}). Default: \code{NULL}.
#' @param assay Character. Name of the assay to run SCTransform on. Default:
#'   \code{"transcript"}.
#' @param verbose Logical. Default: \code{TRUE}.
#' @param ... Additional arguments forwarded to \code{\link[Seurat]{SCTransform}}.
#'
#' @return The updated \code{MatisseObject} with:
#'   \itemize{
#'     \item An \code{"SCT"} assay (created by SCTransform) in the Seurat object.
#'     \item A \code{"pca"} reduction computed on the SCT assay.
#'   }
#'
#' @seealso \code{\link{CreateMatisseObjectFromTranscripts}},
#'   \code{\link[Seurat]{SCTransform}}, \code{\link[Seurat]{RunPCA}}
#'
#' @rdname SCTransformTranscripts
#' @export
setMethod("SCTransformTranscripts", "MatisseObject",
          function(object, n_pca_dims = 50L, vars_to_regress = NULL,
                   assay = "transcript", verbose = TRUE, ...) {
  seu <- object@seurat
  if (is.null(seu[[assay]])) {
    rlang::abort(paste0(
      "No '", assay, "' assay found. ",
      "Run CreateMatisseObjectFromTranscripts() or ",
      "CreateMatisseObject(transcript_counts = ...) first."))
  }

  SeuratObject::DefaultAssay(seu) <- assay

  if (verbose) {
    cli::cli_alert_info(
      "Running SCTransform on '{assay}' assay \\
       ({nrow(seu[[assay]])} features)...")
  }

  seu <- Seurat::SCTransform(
    seu,
    assay           = assay,
    vars.to.regress = vars_to_regress,
    verbose         = verbose,
    ...
  )

  # Cap n_pca_dims to avoid requesting more PCs than available features
  n_pca_dims <- min(as.integer(n_pca_dims),
                    nrow(seu[["SCT"]]) - 1L,
                    ncol(seu) - 1L)

  if (verbose) {
    cli::cli_alert_info("Running PCA with {n_pca_dims} components on SCT assay...")
  }

  seu <- Seurat::RunPCA(seu, assay = "SCT", npcs = n_pca_dims,
                        verbose = verbose)

  object@seurat <- seu

  if (verbose) {
    cli::cli_alert_success(
      "SCTransformTranscripts complete: {n_pca_dims} PCA components computed.")
  }

  object
})
