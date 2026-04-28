#' @include MatisseObject-class.R
#' @include generics.R
NULL

# Manually subset an Assay5 by features. SeuratObject v5's own subset.Assay5
# fails validity ("Layers must be two-dimensional objects") for small feature
# sets and silently reduces layers to vectors in some other paths. This
# helper rebuilds the assay layer-by-layer using .get_assay_layer so the
# v5 vector-reduction bugs are normalised away, then restores meta.features.
.subset_assay_features <- function(assay, features) {
  layer_names <- SeuratObject::Layers(assay)
  if (!"counts" %in% layer_names) {
    rlang::abort("Cannot subset Assay5 features without a 'counts' layer.")
  }
  cnts <- .get_assay_layer(assay, "counts")[features, , drop = FALSE]
  new_a <- SeuratObject::CreateAssay5Object(counts = cnts)
  for (lyr in setdiff(layer_names, "counts")) {
    layer_data <- .get_assay_layer(assay, lyr)
    if (!is.null(layer_data) && length(layer_data) > 0L) {
      new_a <- SeuratObject::SetAssayData(
        new_a, layer = lyr,
        new.data = layer_data[features, , drop = FALSE]
      )
    }
  }
  mf <- assay[[]]
  if (!is.null(mf) && ncol(mf) > 0L) {
    new_a[[]] <- mf[features, , drop = FALSE]
  }
  new_a
}

# SeuratObject v5 bug (verified present in 5.1.0): single-cell subsetting via
# Seurat's [, single_cell] operator reduces Assay5 layers to plain numeric
# vectors with no dim attribute, even though the assay's validity normally
# enforces 2D layers. This helper restores them to sparse matrices using the
# assay's own dim/dimnames so downstream code can index them as matrices.
.get_assay_layer <- function(assay, layer) {
  raw <- methods::slot(assay, "layers")[[layer]]
  if (is.null(raw)) return(NULL)
  if (!is.matrix(raw) && !inherits(raw, "Matrix")) {
    n_features <- nrow(assay)
    n_cells    <- ncol(assay)
    mat <- matrix(raw, nrow = n_features, ncol = n_cells,
                  dimnames = list(rownames(assay), colnames(assay)))
    return(methods::as(mat, "CsparseMatrix"))
  }
  tryCatch(
    SeuratObject::GetAssayData(assay, layer = layer),
    error = function(e) {
      mat <- matrix(as.numeric(raw), nrow = nrow(assay), ncol = ncol(assay),
                    dimnames = list(rownames(assay), colnames(assay)))
      methods::as(mat, "CsparseMatrix")
    }
  )
}

# ---------------------------------------------------------------------------
# show
# ---------------------------------------------------------------------------

#' @describeIn MatisseObject-class Display a summary of a \code{MatisseObject}.
#' @param object A \code{MatisseObject}.
#' @export
setMethod("show", "MatisseObject", function(object) {
  n_cells  <- .n_cells(object)
  n_events <- .n_events(object)

  mode_label <- if (object@input.mode == "junction") "junction-based" else "transcript-based"
  cat("A MatisseObject (", mode_label, " mode)\n", sep = "")
  cat("  Cells        :", n_cells, "\n")
  cat("  Splice events:", n_events, "\n")

  if (!is.null(object@seurat)) {
    n_features <- nrow(object@seurat)
    cat("  Gene features:", n_features, "\n")

    assay_names <- SeuratObject::Assays(object@seurat)
    if (length(assay_names) > 0) {
      default_assay <- SeuratObject::DefaultAssay(object@seurat)
      labelled <- ifelse(assay_names == default_assay,
                         paste0(assay_names, "*"), assay_names)
      cat("  Assays       :", paste(labelled, collapse = ", "), "\n")
    }
    if (length(SeuratObject::Reductions(object@seurat)) > 0) {
      cat("  Reductions   :",
          paste(SeuratObject::Reductions(object@seurat), collapse = ", "), "\n")
    }
  }

  n_isoforms <- .n_isoforms(object)
  if (n_isoforms > 0L) {
    iso_label <- if (object@input.mode == "junction") "Junctions" else "Transcripts"
    cat("  ", formatC(iso_label, width = 13, flag = "-"), ":", n_isoforms, "\n", sep = "")
  }

  # PSI coverage from the "psi" Assay5
  psi_assay <- .get_assay_safe(object@seurat, "psi")
  if (!is.null(psi_assay) && n_events > 0L) {
    psi_ec  <- .get_assay_layer(psi_assay, "data")
    psi_csc <- as(psi_ec, "dgCMatrix")
    n_covered   <- sum(!is.na(psi_csc@x))
    pct_covered <- if (n_cells > 0L && n_events > 0L)
      round(100 * n_covered / (as.double(n_cells) * n_events), 1) else 0
    cat("  PSI coverage :", pct_covered, "% entries covered\n")
  }

  cat("  Version      :", object@version, "\n")
  invisible(object)
})

# ---------------------------------------------------------------------------
# dim / nrow / ncol
# ---------------------------------------------------------------------------

#' @rdname MatisseObject-class
#' @param x A \code{MatisseObject}.
#' @export
setMethod("dim", "MatisseObject", function(x) {
  c(.n_cells(x), .n_events(x))
})

# ---------------------------------------------------------------------------
# Subsetting: [ (cells x events)
# ---------------------------------------------------------------------------

#' @rdname MatisseObject-class
#' @param i Cell barcodes (character) or integer indices.
#' @param j Event IDs (character) or integer indices.
#' @param ... Ignored.
#' @param drop Ignored.
#' @export
setMethod("[", "MatisseObject", function(x, i, j, ..., drop = FALSE) {
  cells_all <- .get_cells(x)

  # Resolve cell indices
  if (missing(i)) {
    cell_names <- cells_all
  } else if (is.character(i)) {
    bad <- setdiff(i, cells_all)
    if (length(bad) > 0) {
      rlang::abort("Some cell barcodes not found in the object.")
    }
    cell_names <- i
  } else {
    cell_names <- cells_all[i]
  }

  # Resolve event indices (only relevant if "psi" assay exists)
  psi_assay <- .get_assay_safe(x@seurat, "psi")
  if (!is.null(psi_assay)) {
    events_all <- rownames(psi_assay)
    if (missing(j)) {
      event_ids <- events_all
    } else if (is.character(j)) {
      bad <- setdiff(j, events_all)
      if (length(bad) > 0) {
        rlang::abort("Some event IDs not found in the PSI assay.")
      }
      event_ids <- j
    } else {
      event_ids <- events_all[j]
    }
  } else {
    event_ids <- character(0)
  }

  # Subset Seurat by cells (handles ALL assays: isoform, psi, gene expression,
  # plus cell metadata and reductions automatically). When event subsetting is
  # requested, also subset the PSI assay's features via .subset_assay_features
  # which normalises around v5's vector-reduction quirks.
  if (!missing(j) && length(event_ids) < 2L) {
    rlang::abort(paste0(
      "Cannot subset a MatisseObject to fewer than 2 events: SeuratObject's ",
      "Assay5 requires at least 2 features. For single-event values, use ",
      "GetPSI(obj)[, event_id] directly."))
  }
  new_seurat <- if (!is.null(x@seurat)) x@seurat[, cell_names] else NULL
  if (!missing(j) && !is.null(new_seurat) && !is.null(psi_assay) &&
      length(event_ids) > 0) {
    new_seurat[["psi"]] <- .subset_assay_features(new_seurat[["psi"]], event_ids)
  }

  methods::new("MatisseObject",
    seurat     = new_seurat,
    input.mode = x@input.mode,
    version    = x@version,
    misc       = x@misc
  )
})

# ---------------------------------------------------------------------------
# Accessor methods
# ---------------------------------------------------------------------------

#' @rdname GetSeurat
setMethod("GetSeurat", "MatisseObject", function(object, ...) object@seurat)

#' @rdname GetPSI
setMethod("GetPSI", "MatisseObject", function(object, ...) {
  psi_assay <- .get_assay_safe(object@seurat, "psi")
  if (is.null(psi_assay)) return(NULL)
  # Seurat: events x cells -> return cells x events (Matisse convention)
  Matrix::t(.get_assay_layer(psi_assay, "data"))
})

#' @rdname SetPSI
setMethod("SetPSI", "MatisseObject", function(object, value) {
  psi_assay <- .get_assay_safe(object@seurat, "psi")
  if (is.null(psi_assay)) {
    rlang::abort("No 'psi' assay exists. Run CalculatePSI() first.")
  }
  # value is cells x events; transpose to events x cells for Seurat
  object@seurat[["psi"]] <- SeuratObject::SetAssayData(
    psi_assay, layer = "data", new.data = Matrix::t(value)
  )
  methods::validObject(object)
  object
})

#' @rdname MatisseMeta
setMethod("MatisseMeta", "MatisseObject", function(object, ...) {
  if (is.null(object@seurat)) return(data.frame())
  object@seurat@meta.data
})

#' @rdname MatisseMeta
setMethod("MatisseMeta<-", "MatisseObject", function(object, value) {
  stopifnot(is.data.frame(value))
  # Delegate to Seurat's AddMetaData so values are aligned by barcode (rownames)
  # rather than by positional index (which the previous loop did, silently
  # writing wrong values when value's rownames disagreed with meta.data's).
  object@seurat <- SeuratObject::AddMetaData(object@seurat, metadata = value)
  object
})

# ---------------------------------------------------------------------------
# [[ operator: metadata first, then Seurat slots
# ---------------------------------------------------------------------------

#' @describeIn MatisseObject-class
#'   Access cell metadata or Seurat slots via \code{[[}.
#'   Returns metadata columns as bare vectors (matching Matisse convention),
#'   then falls back to the embedded Seurat object's \code{[[} for assays,
#'   reductions, etc. The explicit metadata check is load-bearing — Seurat's
#'   own \code{[[} returns metadata columns wrapped in a 1-column data.frame.
#' @aliases [[,MatisseObject-method
#' @export
setMethod("[[", "MatisseObject", function(x, i, j, ...) {
  if (!is.null(x@seurat)) {
    if (i %in% colnames(x@seurat@meta.data)) return(x@seurat@meta.data[[i]])
    return(x@seurat[[i]])
  }
  rlang::abort(paste0("'", i, "' not found in the MatisseObject."))
})

# ---------------------------------------------------------------------------
# $ operator: hybrid dispatch (Seurat metadata -> Seurat/Signac fn -> Seurat [[)
# ---------------------------------------------------------------------------

#' @rdname MatisseObject-class
#' @param name Metadata column name or Seurat/Signac function name.
#' @export
setMethod("$", "MatisseObject", function(x, name) {
  # Priority 1: Seurat cell metadata (most common use: obj$seurat_clusters)
  if (!is.null(x@seurat) && name %in% colnames(x@seurat@meta.data)) {
    return(x@seurat@meta.data[[name]])
  }

  # Priority 2: Seurat or Signac exported function -> return a forwarding closure
  fn <- .find_package_function(name)
  if (!is.null(fn)) {
    force(x)
    force(fn)
    return(function(...) {
      result <- fn(x@seurat, ...)
      if (inherits(result, "Seurat")) {
        x@seurat <- result
        return(x)
      }
      result
    })
  }

  # Priority 3: delegate to Seurat's [[ (handles assays, reductions, etc.)
  if (!is.null(x@seurat)) return(x@seurat[[name]])

  rlang::abort(paste0(
    "'", name, "' not found. Tried (in order): cell metadata column, ",
    "exported function in Seurat or Signac, then Seurat [[ (assays, ",
    "reductions, misc)."))
})

#' @importFrom utils .DollarNames
#' @export
.DollarNames.MatisseObject <- function(x, pattern = "") {
  # Surface both metadata columns AND callable Seurat/Signac functions for
  # tab completion. Without the function names users can't discover that
  # obj$NormalizeData() and similar work via Priority 2 of the $ dispatch.
  meta_nms <- if (!is.null(x@seurat)) colnames(x@seurat@meta.data) else character(0)
  fn_nms <- character(0)
  for (pkg in c("Seurat", "Signac")) {
    pkg_nms <- tryCatch(getNamespaceExports(pkg), error = function(e) character(0))
    fn_nms <- c(fn_nms, pkg_nms)
  }
  grep(pattern, unique(c(meta_nms, fn_nms)), value = TRUE)
}

# Look up an exported function from Seurat or Signac
.find_package_function <- function(name, pkgs = c("Seurat", "Signac")) {
  for (pkg in pkgs) {
    fn <- tryCatch(
      getExportedValue(pkg, name),
      error = function(e) NULL
    )
    if (is.function(fn)) return(fn)
  }
  NULL
}
