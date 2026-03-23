#' Filter a transcript matrix by top-ranked barcodes by retaining the top `cells` barcodes ranked by total counts per barcode.
#'
#' @param tmat A sparse matrix with transcripts/isoforms as rows and cells as columns.
#' @param cells Integer(1). Number of top barcodes with positive counts to retain.
#'
#' @return A filtered sparse matrix containing only the retained barcodes.
#' @export
FilterBarcodes <- function(tmat,
                           cells) {

  # ---- Verify parameter inputs ----
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.", call. = FALSE)
  }

  if (is.null(tmat) || !(inherits(tmat, "Matrix") || is.matrix(tmat))) {
    stop("tmat must be a matrix-like object with cells as columns.", call. = FALSE)
  }

  if (length(cells) != 1L || is.na(cells) || cells < 1L || cells != as.integer(cells)) {
    stop("cells must be a positive integer(1).", call. = FALSE)
  }
  cells <- as.integer(cells)

  if (!inherits(tmat, "dgCMatrix")) {
    tmat <- methods::as(tmat, "dgCMatrix")
  }

  total.counts <- Matrix::colSums(tmat)
  keep.nz <- total.counts > 0

  if (!any(keep.nz)) {
    stop("No barcodes with positive total counts were found.", call. = FALSE)
  }

  tmat <- tmat[, keep.nz, drop = FALSE]
  total.counts <- total.counts[keep.nz]

  ord <- order(total.counts, decreasing = TRUE)

  if (cells > length(ord)) {
    stop("cells exceeds the number of barcodes with positive counts.", call. = FALSE)
  }

  # ---- Filter transcript matrix ----
  keep <- ord[seq_len(cells)]
  out <- tmat[, keep, drop = FALSE]

  attr(out, "barcode.cutoff") <- unname(sort(total.counts, decreasing = TRUE)[cells])
  attr(out, "retained.barcodes") <- colnames(out)

  out
}
