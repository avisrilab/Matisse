#' @include MatisseObject-class.R
#' @include generics.R
NULL

# ---------------------------------------------------------------------------
# show
# ---------------------------------------------------------------------------

#' @describeIn MatisseObject-class Display a summary of a \code{MatisseObject}.
#' @param object A \code{MatisseObject}.
#' @export
setMethod("show", "MatisseObject", function(object) {
  n_cells  <- .n_cells(object)
  n_events <- .n_events(object)

  cat("A MatisseObject\n")
  cat("  Cells        :", n_cells, "\n")
  cat("  Splice events:", n_events, "\n")

  # Seurat summary line
  if (!is.null(object@seurat)) {
    n_features <- nrow(object@seurat)
    cat("  Gene features:", n_features, "\n")
    if (length(SeuratObject::Reductions(object@seurat)) > 0) {
      cat("  Reductions   :",
          paste(SeuratObject::Reductions(object@seurat), collapse = ", "), "\n")
    }
  }

  # Junction layer
  if (!is.null(object@junction_counts)) {
    cat("  Junctions    :", ncol(object@junction_counts), "\n")
  }

  # PSI coverage — count non-NA stored entries without dense coercion
  if (!is.null(object@psi)) {
    psi_csc     <- as(object@psi, "dgCMatrix")
    n_covered   <- sum(!is.na(psi_csc@x))
    pct_covered <- if (n_cells > 0L && n_events > 0L)
      round(100 * n_covered / (n_cells * n_events), 1) else 0
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
  # Resolve cell indices
  cells_all <- .get_cells(x)
  if (missing(i)) {
    cell_idx <- seq_along(cells_all)
  } else if (is.character(i)) {
    cell_idx <- match(i, cells_all)
    if (any(is.na(cell_idx))) {
      rlang::abort("Some cell barcodes not found in the object.")
    }
  } else {
    cell_idx <- i
  }

  # Resolve event indices
  if (!is.null(x@psi)) {
    events_all <- colnames(x@psi)
    if (missing(j)) {
      event_idx <- seq_along(events_all)
    } else if (is.character(j)) {
      event_idx <- match(j, events_all)
      if (any(is.na(event_idx))) {
        rlang::abort("Some event IDs not found in the PSI matrix.")
      }
    } else {
      event_idx <- j
    }
  } else {
    event_idx <- integer(0)
  }

  # Subset Seurat object
  new_seurat <- if (!is.null(x@seurat)) {
    x@seurat[, cell_idx]
  } else {
    NULL
  }

  # Subset matrices
  subset_mat <- function(m) {
    if (is.null(m)) return(NULL)
    m[cell_idx, event_idx, drop = FALSE]
  }
  subset_jxn <- function(m) {
    if (is.null(m)) return(NULL)
    m[cell_idx, , drop = FALSE]
  }

  # Subset isoform metadata
  new_meta <- if (nrow(x@isoform_metadata) > 0) {
    x@isoform_metadata[cell_idx, , drop = FALSE]
  } else {
    data.frame()
  }

  # Subset event_data
  new_event_data <- if (nrow(x@event_data) > 0 && length(event_idx) > 0) {
    x@event_data[event_idx, , drop = FALSE]
  } else {
    x@event_data
  }

  methods::new("MatisseObject",
    seurat           = new_seurat,
    psi              = subset_mat(x@psi),
    inclusion_counts = subset_mat(x@inclusion_counts),
    exclusion_counts = subset_mat(x@exclusion_counts),
    junction_counts  = subset_jxn(x@junction_counts),
    event_data       = new_event_data,
    junction_data    = x@junction_data,
    isoform_metadata = new_meta,
    version          = x@version,
    misc             = x@misc
  )
})

# ---------------------------------------------------------------------------
# Accessor methods
# ---------------------------------------------------------------------------

#' @rdname GetSeurat
setMethod("GetSeurat", "MatisseObject", function(object, ...) object@seurat)
#' @rdname GetPSI
setMethod("GetPSI",    "MatisseObject", function(object, ...) object@psi)

#' @rdname SetPSI
setMethod("SetPSI", "MatisseObject", function(object, value) {
  object@psi <- value
  methods::validObject(object)
  object
})

#' @rdname GetJunctionCounts
setMethod("GetJunctionCounts",  "MatisseObject",
          function(object, ...) object@junction_counts)
#' @rdname GetInclusionCounts
setMethod("GetInclusionCounts", "MatisseObject",
          function(object, ...) object@inclusion_counts)
#' @rdname GetExclusionCounts
setMethod("GetExclusionCounts", "MatisseObject",
          function(object, ...) object@exclusion_counts)
#' @rdname GetEventData
setMethod("GetEventData",       "MatisseObject",
          function(object, ...) object@event_data)
#' @rdname GetJunctionData
setMethod("GetJunctionData",    "MatisseObject",
          function(object, ...) object@junction_data)
#' @rdname MatisseMeta
setMethod("MatisseMeta",        "MatisseObject",
          function(object, ...) object@isoform_metadata)

#' @rdname MatisseMeta
setMethod("MatisseMeta<-", "MatisseObject", function(object, value) {
  stopifnot(is.data.frame(value))
  object@isoform_metadata <- value
  object
})

#' @describeIn AddIsoformMetadata Add or update columns in the isoform metadata.
setMethod("AddIsoformMetadata", "MatisseObject",
          function(object, metadata, ...) {
  cells <- .get_cells(object)

  if (is.data.frame(metadata)) {
    new_cols <- metadata
  } else if (is.numeric(metadata) || is.character(metadata) ||
             is.logical(metadata) || is.integer(metadata)) {
    if (is.null(names(metadata))) {
      rlang::abort("'metadata' vector must be named with cell barcodes.")
    }
    # as.data.frame() on a named vector produces N-row x 1-col with row names
    new_cols <- as.data.frame(metadata)
  } else {
    rlang::abort("'metadata' must be a data.frame or named vector.")
  }

  # Align rows by barcode
  if (!is.null(rownames(new_cols))) {
    new_cols <- new_cols[cells, , drop = FALSE]
  }

  current <- object@isoform_metadata
  if (nrow(current) == 0) {
    current <- data.frame(row.names = cells)
  }

  for (col in colnames(new_cols)) {
    current[[col]] <- new_cols[[col]]
  }

  object@isoform_metadata <- current
  object
})

# ---------------------------------------------------------------------------
# [[ operator: isoform metadata first, then Seurat
# ---------------------------------------------------------------------------

#' @describeIn MatisseObject-class
#'   Access isoform metadata columns or Seurat slots via \code{[[}.
#'   Checks \code{isoform_metadata} first; falls back to the embedded
#'   Seurat object.
#' @aliases [[,MatisseObject-method
#' @export
setMethod("[[", "MatisseObject", function(x, i, j, ...) {
  meta <- x@isoform_metadata
  if (nrow(meta) > 0 && i %in% colnames(meta)) {
    return(meta[[i]])
  }
  if (!is.null(x@seurat)) {
    return(x@seurat[[i]])
  }
  rlang::abort(paste0(
    "'", i, "' not found in isoform_metadata or the embedded Seurat object."))
})

# ---------------------------------------------------------------------------
# $ operator: forwards to Seurat
# ---------------------------------------------------------------------------

#' @rdname MatisseObject-class
#' @param name Column name in the Seurat \code{meta.data} (used by \code{$}).
#' @export
setMethod("$", "MatisseObject", function(x, name) {
  if (is.null(x@seurat)) {
    rlang::abort("No Seurat object is embedded in this MatisseObject.")
  }
  x@seurat[[name]]
})
