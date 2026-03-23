#' Resolve bin-level support to transcript-level support.
#'
#' This function resolves both bin-level signal and bin-level read support to the
#' transcript level using transcript quantification weights for shared bins.
#'
#' @param model List produced by `BuildGeneModel()`. Must contain `bins`,
#'   `transcript_bin_map`, and `transcript_table`.
#' @param obs List produced by `ComputeBinCoverage()`. Must contain
#'   `signal_matrix` and `read_matrix`.
#' @param object A Seurat object containing transcript quantifications.
#' @param transcript.assay Character(1). Assay name for transcript quantifications.
#' @param quant.layer Character(1). Layer name to retrieve from the assay.
#' @param group.by Optional Character(1). Metadata column in `object` to append to
#'   the returned support tables.
#'
#' @return A list containing unique-only and resolved transcript-level signal and
#'   read-support tables.
#' @export
ResolveTranscriptCoverage <- function(
    model,
    obs,
    object,
    transcript.assay = "transcript",
    quant.layer = "data",
    group.by = NULL
) {

  # ---- Helper checks ----
  .check_named_list <- function(x, nm) {
    if (!is.list(x)) {
      stop(nm, " must be a list.", call. = FALSE)
    }
  }

  .check_required_fields <- function(x, fields, nm) {
    miss <- setdiff(fields, names(x))
    if (length(miss) > 0L) {
      stop(
        nm, " is missing required field(s): ",
        paste(miss, collapse = ", "),
        call. = FALSE
      )
    }
  }

  .check_scalar_character <- function(x, nm) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
      stop(nm, " must be a non-empty character(1).", call. = FALSE)
    }
  }

  .check_matrix_like <- function(x, nm) {
    if (!(methods::is(x, "matrix") || methods::is(x, "Matrix"))) {
      stop(nm, " must be a matrix-like object.", call. = FALSE)
    }
    if (is.null(rownames(x)) || anyNA(rownames(x)) || any(!nzchar(rownames(x)))) {
      stop(nm, " must have non-empty rownames corresponding to cell IDs.", call. = FALSE)
    }
    if (is.null(colnames(x)) || anyNA(colnames(x)) || any(!nzchar(colnames(x)))) {
      stop(nm, " must have non-empty colnames corresponding to bin IDs.", call. = FALSE)
    }
  }

  .bind_or_empty <- function(x) {
    keep <- !vapply(x, is.null, logical(1))
    if (!any(keep)) {
      return(data.frame(
        cell = character(0),
        transcript_id = character(0),
        bin_id = character(0),
        value = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    out <- do.call(rbind, x[keep])
    rownames(out) <- NULL
    out
  }

  # ---- Validate inputs ----
  .check_named_list(model, "model")
  .check_required_fields(
    model,
    c("bins", "transcript_bin_map", "transcript_table"),
    "model"
  )

  .check_named_list(obs, "obs")
  .check_required_fields(obs, c("signal_matrix", "read_matrix"), "obs")

  .check_scalar_character(transcript.assay, "transcript.assay")
  .check_scalar_character(quant.layer, "quant.layer")

  if (!is.null(group.by)) {
    .check_scalar_character(group.by, "group.by")
  }

  if (!methods::is(model$bins, "GRanges")) {
    stop("model$bins must be a GenomicRanges::GRanges object.", call. = FALSE)
  }

  if (!is.data.frame(model$transcript_bin_map)) {
    stop("model$transcript_bin_map must be a data.frame.", call. = FALSE)
  }

  if (!all(c("tx_id", "bin_id") %in% colnames(model$transcript_bin_map))) {
    stop(
      "model$transcript_bin_map must contain columns 'tx_id' and 'bin_id'.",
      call. = FALSE
    )
  }

  if (!is.data.frame(model$transcript_table)) {
    stop("model$transcript_table must be a data.frame.", call. = FALSE)
  }

  if (!all(c("tx_id", "tx_name") %in% colnames(model$transcript_table))) {
    stop(
      "model$transcript_table must contain columns 'tx_id' and 'tx_name'.",
      call. = FALSE
    )
  }

  signal.mat <- obs$signal_matrix
  read.mat <- obs$read_matrix

  .check_matrix_like(signal.mat, "obs$signal_matrix")
  .check_matrix_like(read.mat, "obs$read_matrix")

  if (!identical(rownames(signal.mat), rownames(read.mat))) {
    stop("obs$signal_matrix and obs$read_matrix must have identical rownames.", call. = FALSE)
  }

  if (!identical(colnames(signal.mat), colnames(read.mat))) {
    stop("obs$signal_matrix and obs$read_matrix must have identical colnames.", call. = FALSE)
  }

  bin.mcols <- S4Vectors::mcols(model$bins)
  if (!("bin_id" %in% colnames(bin.mcols))) {
    stop("model$bins must contain a 'bin_id' metadata column.", call. = FALSE)
  }

  if (anyNA(bin.mcols$bin_id) || any(!nzchar(as.character(bin.mcols$bin_id)))) {
    stop("model$bins$bin_id contains missing or empty values.", call. = FALSE)
  }

  if (anyNA(model$transcript_bin_map$tx_id) || any(!nzchar(as.character(model$transcript_bin_map$tx_id)))) {
    stop("model$transcript_bin_map$tx_id contains missing or empty values.", call. = FALSE)
  }

  if (anyNA(model$transcript_bin_map$bin_id) || any(!nzchar(as.character(model$transcript_bin_map$bin_id)))) {
    stop("model$transcript_bin_map$bin_id contains missing or empty values.", call. = FALSE)
  }

  if (anyNA(model$transcript_table$tx_id) || any(!nzchar(as.character(model$transcript_table$tx_id)))) {
    stop("model$transcript_table$tx_id contains missing or empty values.", call. = FALSE)
  }

  if (anyNA(model$transcript_table$tx_name) || any(!nzchar(as.character(model$transcript_table$tx_name)))) {
    stop("model$transcript_table$tx_name contains missing or empty values.", call. = FALSE)
  }

  if (!isS4(object)) {
    stop("object must be a Seurat object.", call. = FALSE)
  }

  assay.names <- names(object@assays)
  if (!(transcript.assay %in% assay.names)) {
    stop("transcript.assay not found in object: ", transcript.assay, call. = FALSE)
  }

  if (!is.null(group.by)) {
    object.meta <- object[[]]
    if (!(group.by %in% colnames(object.meta))) {
      stop("group.by not found in object metadata: ", group.by, call. = FALSE)
    }
  }

  # ---- Access inputs ----
  bins <- model$bins
  tx.bin.map <- model$transcript_bin_map
  tx.df <- model$transcript_table

  # ---- Access transcript quantifications ----
  quant.mat <- Seurat::GetAssayData(
    object = object,
    assay = transcript.assay,
    layer = quant.layer
  )
  quant.mat <- as.matrix(quant.mat)

  if (is.null(rownames(quant.mat)) || anyNA(rownames(quant.mat)) || any(!nzchar(rownames(quant.mat)))) {
    stop("Transcript quantification matrix must have non-empty rownames corresponding to transcript IDs.", call. = FALSE)
  }

  if (is.null(colnames(quant.mat)) || anyNA(colnames(quant.mat)) || any(!nzchar(colnames(quant.mat)))) {
    stop("Transcript quantification matrix must have non-empty colnames corresponding to cell IDs.", call. = FALSE)
  }

  # ---- Restrict to shared cells ----
  common.cells <- intersect(rownames(signal.mat), colnames(quant.mat))
  signal.mat <- signal.mat[common.cells, , drop = FALSE]
  read.mat <- read.mat[common.cells, , drop = FALSE]
  quant.mat <- quant.mat[, common.cells, drop = FALSE]

  # ---- Restrict transcript metadata/map to transcripts present in quantification matrix ----
  tx.df <- tx.df[tx.df$tx_id %in% rownames(quant.mat), , drop = FALSE]
  tx.bin.map <- tx.bin.map[tx.bin.map$tx_id %in% rownames(quant.mat), , drop = FALSE]

  bin.ids <- colnames(signal.mat)
  cell.ids <- rownames(signal.mat)

  # ---- Build bin -> transcripts lookup ----
  bin.to.tx <- split(tx.bin.map$tx_id, tx.bin.map$bin_id)
  bin.to.tx <- lapply(bin.to.tx, unique)

  # ---- Build output containers ----
  unique.signal.list <- vector("list", length(bin.ids))
  resolved.signal.list <- vector("list", length(bin.ids))
  unique.read.list <- vector("list", length(bin.ids))
  resolved.read.list <- vector("list", length(bin.ids))

  names(unique.signal.list) <- bin.ids
  names(resolved.signal.list) <- bin.ids
  names(unique.read.list) <- bin.ids
  names(resolved.read.list) <- bin.ids

  # ---- Resolve each bin ----
  for (bin.id in bin.ids) {
    txs <- bin.to.tx[[bin.id]]

    if (length(txs) == 0L) {
      next
    }

    y.signal <- as.numeric(signal.mat[, bin.id])
    y.read <- as.numeric(read.mat[, bin.id])
    names(y.signal) <- cell.ids
    names(y.read) <- cell.ids

    if (all(y.signal == 0) && all(y.read == 0)) {
      next
    }

    # ---- Unique-only assignment ----
    if (length(txs) == 1L) {
      unique.signal.list[[bin.id]] <- data.frame(
        cell = cell.ids,
        transcript_id = txs,
        bin_id = bin.id,
        value = y.signal,
        stringsAsFactors = FALSE
      )

      resolved.signal.list[[bin.id]] <- data.frame(
        cell = cell.ids,
        transcript_id = txs,
        bin_id = bin.id,
        value = y.signal,
        stringsAsFactors = FALSE
      )

      unique.read.list[[bin.id]] <- data.frame(
        cell = cell.ids,
        transcript_id = txs,
        bin_id = bin.id,
        value = y.read,
        stringsAsFactors = FALSE
      )

      resolved.read.list[[bin.id]] <- data.frame(
        cell = cell.ids,
        transcript_id = txs,
        bin_id = bin.id,
        value = y.read,
        stringsAsFactors = FALSE
      )

    } else {
      # ---- Quantification-resolved assignment for shared bins ----
      txs.use <- txs[txs %in% rownames(quant.mat)]

      if (length(txs.use) == 0L) {
        next
      }

      q <- quant.mat[txs.use, cell.ids, drop = FALSE]
      denom <- colSums(q)

      weight.mat <- matrix(
        0,
        nrow = nrow(q),
        ncol = ncol(q),
        dimnames = dimnames(q)
      )

      nonzero <- denom > 0
      if (any(nonzero)) {
        weight.mat[, nonzero] <- sweep(
          q[, nonzero, drop = FALSE],
          2,
          denom[nonzero],
          "/"
        )
      }

      signal.value.mat <- sweep(weight.mat, 2, y.signal, "*")
      read.value.mat <- sweep(weight.mat, 2, y.read, "*")

      resolved.signal.df <- as.data.frame(as.table(signal.value.mat), stringsAsFactors = FALSE)
      colnames(resolved.signal.df) <- c("transcript_id", "cell", "value")
      resolved.signal.df$bin_id <- bin.id
      resolved.signal.df <- resolved.signal.df[, c("cell", "transcript_id", "bin_id", "value")]
      resolved.signal.list[[bin.id]] <- resolved.signal.df

      resolved.read.df <- as.data.frame(as.table(read.value.mat), stringsAsFactors = FALSE)
      colnames(resolved.read.df) <- c("transcript_id", "cell", "value")
      resolved.read.df$bin_id <- bin.id
      resolved.read.df <- resolved.read.df[, c("cell", "transcript_id", "bin_id", "value")]
      resolved.read.list[[bin.id]] <- resolved.read.df
    }
  }

  # ---- Combine outputs ----
  unique.signal <- .bind_or_empty(unique.signal.list)
  resolved.signal <- .bind_or_empty(resolved.signal.list)
  unique.read <- .bind_or_empty(unique.read.list)
  resolved.read <- .bind_or_empty(resolved.read.list)

  # ---- Add transcript names ----
  tx.name.map <- tx.df$tx_name
  names(tx.name.map) <- tx.df$tx_id

  unique.signal$transcript_name <- tx.name.map[unique.signal$transcript_id]
  resolved.signal$transcript_name <- tx.name.map[resolved.signal$transcript_id]
  unique.read$transcript_name <- tx.name.map[unique.read$transcript_id]
  resolved.read$transcript_name <- tx.name.map[resolved.read$transcript_id]

  # ---- Add grouping metadata if requested ----
  if (!is.null(group.by)) {
    group.vec <- object[[group.by]][, 1]
    names(group.vec) <- colnames(object)

    unique.signal$group <- group.vec[unique.signal$cell]
    resolved.signal$group <- group.vec[resolved.signal$cell]
    unique.read$group <- group.vec[unique.read$cell]
    resolved.read$group <- group.vec[resolved.read$cell]
  }

  # ---- Return ----
  list(
    unique_signal = unique.signal,
    resolved_signal = resolved.signal,
    unique_read_support = unique.read,
    resolved_read_support = resolved.read,
    cells = common.cells,
    bins = bins,
    transcript_table = tx.df
  )
}
