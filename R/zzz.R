# ---------------------------------------------------------------------------
# Package load hooks
# ---------------------------------------------------------------------------
# At load time we copy the formals of each underlying Seurat / Signac S3
# method onto the matching `Foo.MatisseObject` dispatcher. This gives
# callers proper tab completion and argument hints when they write
# `FindMarkers(matisse_obj, |)` etc., without us hand-mirroring every
# Seurat signature in dispatch.R. The body installed alongside the new
# formals captures all bound formals plus `...`, swaps `object` for the
# embedded Seurat object, and forwards via do.call -- so named args
# defined as part of the copied signature reach the underlying function
# correctly (a plain `...`-only body would drop them).
#
# Dispatchers with custom logic that shouldn't be replaced (e.g.
# SCTransform.MatisseObject's mode-aware default assay) are intentionally
# omitted from the jobs list below.

.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)

  make_body <- function(seurat_fn_expr) {
    bquote({
      .args        <- mget(
        setdiff(names(formals(sys.function())), "..."),
        envir = environment()
      )
      .matisse_obj <- .args$object
      .args$object <- .matisse_obj@seurat
      .result      <- do.call(.(seurat_fn_expr), c(.args, list(...)))
      if (inherits(.result, "Seurat")) {
        .matisse_obj@seurat <- .result
        return(.matisse_obj)
      }
      .result
    })
  }

  jobs <- list(
    list(name = "RunPCA",          pkg = "Seurat",       class = "Seurat",
         fn = quote(Seurat::RunPCA)),
    list(name = "RunUMAP",         pkg = "Seurat",       class = "Seurat",
         fn = quote(Seurat::RunUMAP)),
    list(name = "FindNeighbors",   pkg = "Seurat",       class = "Seurat",
         fn = quote(Seurat::FindNeighbors)),
    list(name = "FindClusters",    pkg = "Seurat",       class = "Seurat",
         fn = quote(Seurat::FindClusters)),
    list(name = "FindMarkers",     pkg = "Seurat",       class = "Seurat",
         fn = quote(Seurat::FindMarkers)),
    list(name = "AddMetaData",     pkg = "SeuratObject", class = "Seurat",
         fn = quote(SeuratObject::AddMetaData)),
    list(name = "RunTFIDF",        pkg = "Signac",       class = "Assay",
         fn = quote(Signac::RunTFIDF)),
    list(name = "RunSVD",          pkg = "Signac",       class = "Assay",
         fn = quote(Signac::RunSVD)),
    list(name = "FindTopFeatures", pkg = "Signac",       class = "Assay",
         fn = quote(Signac::FindTopFeatures))
  )

  for (j in jobs) {
    target_name <- paste0(j$name, ".MatisseObject")
    if (!exists(target_name, envir = ns, inherits = FALSE)) next
    src <- tryCatch(
      utils::getS3method(j$name, j$class, optional = TRUE,
                          envir = asNamespace(j$pkg)),
      error = function(e) NULL
    )
    if (is.null(src) || !is.function(src)) next
    target       <- get(target_name, envir = ns)
    src_formals  <- formals(src)
    # The body uses `list(...)` to forward extra args, which requires the
    # function to declare `...`. Most Seurat S3 methods do, but
    # SeuratObject::AddMetaData.Seurat is `(object, metadata, col.name)`
    # with no dots. Append `...` if it's missing.
    if (!"..." %in% names(src_formals)) {
      src_formals <- c(src_formals, alist(... = ))
    }
    formals(target) <- src_formals
    body(target)    <- make_body(j$fn)
    environment(target) <- ns
    assignInNamespace(target_name, target, ns = pkgname)
  }
}
