# ---------------------------------------------------------------------------
# Package load hooks
# ---------------------------------------------------------------------------
# At load time we copy the formals of each underlying Seurat / Signac S3
# method onto the matching `Foo.MatisseObject` dispatcher. This gives
# callers proper tab completion and argument hints when they write
# `FindMarkers(matisse_obj, |)` etc., without us hand-mirroring every
# Seurat signature in dispatch.R.
#
# The body installed alongside the new formals is a *literal call* that
# forwards each named formal by name, plus `...`. Critically it does NOT
# pack the embedded Seurat object into an intermediate list (via
# `mget` / `do.call`, or `match.call` / `eval`). Doing so bumps R's NAMED
# reference count on the Seurat, which causes Seurat's internal slot
# mutations (e.g. `object@reductions <- ...`) to deep-copy the entire
# object on every assignment. Empirically that path was 8-15x slower than
# a direct call on RunPCA / RunUMAP at typical single-cell scales.
#
# Dispatchers with custom logic that shouldn't be replaced (e.g.
# SCTransform.MatisseObject's mode-aware default assay) are intentionally
# omitted from the jobs list below.

.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)

  # Build a body of the form:
  #   {
  #     .matisse_obj <- object
  #     .result <- <fn>(object = object@seurat, fmA = fmA, fmB = fmB, ...)
  #     if (inherits(.result, "Seurat")) { .matisse_obj@seurat <- .result; return(.matisse_obj) }
  #     .result
  #   }
  # Every named formal of the wrapped Seurat / Signac method is forwarded
  # by symbol reference. R's missing-propagation handles formals the user
  # didn't pass that also have no default (e.g. AddMetaData's `metadata`):
  # `f(x = x)` with a missing `x` passes the missingness through, so the
  # receiving function errors with its own "argument 'x' is missing"
  # message rather than silently dropping the arg.
  make_body <- function(seurat_fn_expr, src_formals) {
    arg_names_all <- names(src_formals)
    arg_names     <- setdiff(arg_names_all, c("object", "..."))
    has_dots      <- "..." %in% arg_names_all
    arg_exprs <- c(
      list(object = quote(object@seurat)),
      stats::setNames(lapply(arg_names, as.name), arg_names),
      if (has_dots) list(quote(...)) else list()
    )
    call_expr <- as.call(c(list(seurat_fn_expr), arg_exprs))
    bquote({
      .matisse_obj <- object
      .result      <- .(call_expr)
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
    # Most Seurat S3 methods declare `...`, but SeuratObject::AddMetaData.Seurat
    # is `(object, metadata, col.name)` with no dots. Append `...` if missing
    # so users can still pass through extra named args.
    if (!"..." %in% names(src_formals)) {
      src_formals <- c(src_formals, alist(... = ))
    }
    formals(target) <- src_formals
    body(target)    <- make_body(j$fn, src_formals)
    environment(target) <- ns
    assignInNamespace(target_name, target, ns = pkgname)
  }
}
