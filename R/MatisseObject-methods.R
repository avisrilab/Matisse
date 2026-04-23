#' @include MatisseObject-class.R
#' @include generics.R
NULL

# SeuratObject v5 bug: single-cell subsetting reduces Assay5 layers to vectors.
# This helper converts them back to sparse matrices.
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

  mode_label <- if (object@mode == "junction") "junction-based" else "event-based"
  cat("A MatisseObject (", mode_label, " mode)\n", sep = "")
  cat("  Cells        :", n_cells, "\n")
  cat("  Splice events:", n_events, "\n")

  if (!is.null(object@seurat)) {
    n_features <- nrow(object@seurat)
    cat("  Gene features:", n_features, "\n")

    assay_names <- SeuratObject::Assays(object@seurat)
    if (length(assay_names) > 0) {
      cat("  Assays       :", paste(assay_names, collapse = ", "), "\n")
    }
    if (length(SeuratObject::Reductions(object@seurat)) > 0) {
      cat("  Reductions   :",
          paste(SeuratObject::Reductions(object@seurat), collapse = ", "), "\n")
    }
  }

  n_junctions <- .n_junctions(object)
  if (n_junctions > 0L) {
    cat("  Junctions    :", n_junctions, "\n")
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

  # Subset Seurat by cells (handles ALL assays: junction, psi, transcript,
  # gene expression, plus cell metadata and reductions automatically)
  new_seurat <- if (!is.null(x@seurat)) x@seurat[, cell_names] else NULL

  # Subset event_data to selected events (PSI assay keeps all features;
  # GetPSI/GetInclusionCounts/GetExclusionCounts filter by event_data)
  if (missing(j)) {
    new_event_data <- x@event_data
  } else if (length(event_ids) == 0L) {
    new_event_data <- x@event_data[integer(0), , drop = FALSE]
  } else if (nrow(x@event_data) > 0) {
    new_event_data <- x@event_data[match(event_ids, x@event_data$event_id), , drop = FALSE]
  } else {
    new_event_data <- x@event_data
  }

  methods::new("MatisseObject",
    seurat        = new_seurat,
    event_data    = new_event_data,
    junction_data = x@junction_data,
    mode          = x@mode,
    version       = x@version,
    misc          = x@misc
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
  # Seurat: events x cells → return cells x events (Matisse convention)
  psi_ec  <- .get_assay_layer(psi_assay, "data")
  psi_ce  <- Matrix::t(psi_ec)
  # Filter to events active in event_data (assay may hold a superset)
  if (nrow(object@event_data) > 0) {
    active <- intersect(object@event_data$event_id, colnames(psi_ce))
    psi_ce <- psi_ce[, active, drop = FALSE]
  }
  psi_ce
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

#' @rdname GetJunctionCounts
setMethod("GetJunctionCounts", "MatisseObject", function(object, ...) {
  jxn_assay <- .get_assay_safe(object@seurat, "junction")
  if (is.null(jxn_assay)) return(NULL)
  # Assay5 is junctions x cells; return cells x junctions (Matisse convention)
  jxn_ec <- .get_assay_layer(jxn_assay, "counts")
  Matrix::t(jxn_ec)
})

#' @rdname GetInclusionCounts
setMethod("GetInclusionCounts", "MatisseObject", function(object, ...) {
  psi_assay <- .get_assay_safe(object@seurat, "psi")
  if (is.null(psi_assay)) return(NULL)
  inc_ec <- .get_assay_layer(psi_assay, "counts")
  inc_ce <- Matrix::t(inc_ec)
  if (nrow(object@event_data) > 0) {
    active <- intersect(object@event_data$event_id, colnames(inc_ce))
    inc_ce <- inc_ce[, active, drop = FALSE]
  }
  inc_ce
})

#' @rdname GetExclusionCounts
setMethod("GetExclusionCounts", "MatisseObject", function(object, ...) {
  psi_assay <- .get_assay_safe(object@seurat, "psi")
  if (is.null(psi_assay)) return(NULL)
  exc_ec <- .get_assay_layer(psi_assay, "exclusion")
  if (is.null(exc_ec) || length(exc_ec) == 0) return(NULL)
  exc_ce <- Matrix::t(exc_ec)
  if (nrow(object@event_data) > 0) {
    active <- intersect(object@event_data$event_id, colnames(exc_ce))
    exc_ce <- exc_ce[, active, drop = FALSE]
  }
  exc_ce
})

#' @rdname GetTranscriptCounts
setMethod("GetTranscriptCounts", "MatisseObject", function(object, ...) {
  tx_assay <- .get_assay_safe(object@seurat, "transcript")
  if (is.null(tx_assay)) return(NULL)
  # Already in transcripts x cells (Seurat convention); return as-is
  SeuratObject::GetAssayData(tx_assay, layer = "counts")
})

#' @rdname GetEventData
setMethod("GetEventData", "MatisseObject",
          function(object, ...) object@event_data)

#' @rdname GetJunctionData
setMethod("GetJunctionData", "MatisseObject",
          function(object, ...) object@junction_data)

#' @rdname MatisseMeta
setMethod("MatisseMeta", "MatisseObject", function(object, ...) {
  if (is.null(object@seurat)) return(data.frame())
  object@seurat@meta.data
})

#' @rdname MatisseMeta
setMethod("MatisseMeta<-", "MatisseObject", function(object, value) {
  stopifnot(is.data.frame(value))
  # Merge new columns into seurat@meta.data rather than replacing everything
  for (col in colnames(value)) {
    object@seurat@meta.data[[col]] <- value[[col]]
  }
  object
})

#' @describeIn AddIsoformMetadata Add or update columns in the cell metadata.
setMethod("AddIsoformMetadata", "MatisseObject",
          function(object, metadata, ...) {
  # Delegate to Seurat's AddMetaData — it handles data.frames and named vectors
  result <- SeuratObject::AddMetaData(object@seurat, metadata = metadata)
  if (inherits(result, "Seurat")) object@seurat <- result
  object
})

# ---------------------------------------------------------------------------
# [[ operator: metadata first, then Seurat slots
# ---------------------------------------------------------------------------

#' @describeIn MatisseObject-class
#'   Access cell metadata or Seurat slots via \code{[[}.
#'   Checks \code{seurat@@meta.data} first; falls back to the embedded
#'   Seurat object (assays, reductions, etc.).
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
# $ operator: hybrid dispatch (Seurat metadata → Seurat/Signac fn → Seurat [[)
# ---------------------------------------------------------------------------

#' @rdname MatisseObject-class
#' @param name Metadata column name or Seurat/Signac function name.
#' @export
setMethod("$", "MatisseObject", function(x, name) {
  # Priority 1: Seurat cell metadata (most common use: obj$seurat_clusters)
  if (!is.null(x@seurat) && name %in% colnames(x@seurat@meta.data)) {
    return(x@seurat@meta.data[[name]])
  }

  # Priority 2: Seurat or Signac exported function → return a forwarding closure
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

  rlang::abort(paste0("'", name, "' not found in Seurat metadata, ",
    "Seurat/Signac functions, or the Seurat object."))
})

#' @export
.DollarNames.MatisseObject <- function(x, pattern = "") {
  nms <- if (!is.null(x@seurat)) colnames(x@seurat@meta.data) else character(0)
  grep(pattern, nms, value = TRUE)
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
