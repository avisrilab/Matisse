#' Visualize a barcode-rank plot from a transcript count matrix
#'
#' @param tmat A count matrix with transcripts/isoforms as rows and cells as columns.
#' @param cells Integer(1). Number of top-ranked barcodes to highlight.
#' @param verbose Logical(1). Whether to print the selected rank and cutoff.
#'
#' @return Invisibly returns a data.frame with columns `barcode`, `total.counts`, and `rank`.
#' @export
PlotBarcodeRank <- function(tmat, 
                            cells, 
                            verbose = TRUE) {
  
  # ---- Verify parameter inputs ----
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.", call. = FALSE)
  }
  
  if (is.null(tmat) || !(inherits(tmat, "Matrix") || is.matrix(tmat))) {
    stop("tmat must be a matrix-like object with cells in columns.", call. = FALSE)
  }
  
  if (is.null(colnames(tmat)) || anyNA(colnames(tmat)) || any(!nzchar(colnames(tmat)))) {
    stop("tmat must have non-empty colnames containing barcodes.", call. = FALSE)
  }
  
  if (length(cells) != 1L || is.na(cells) || cells < 1L || cells != as.integer(cells)) {
    stop("cells must be a positive integer(1).", call. = FALSE)
  }
  cells <- as.integer(cells)
  
  if (length(verbose) != 1L || is.na(verbose) || !is.logical(verbose)) {
    stop("verbose must be a logical(1).", call. = FALSE)
  }
  
  if (!inherits(tmat, "dgCMatrix")) {
    tmat <- methods::as(tmat, "dgCMatrix")
  }
  
  # ---- Compute counts ----
  total.counts <- Matrix::colSums(tmat)
  barcodes <- as.character(colnames(tmat))
  
  keep <- total.counts > 0
  if (!any(keep)) {
    stop("No barcodes with positive total counts were found.", call. = FALSE)
  }
  
  total.counts <- total.counts[keep]
  barcodes <- barcodes[keep]
  
  ord <- order(total.counts, decreasing = TRUE)
  
  barcode.df <- data.frame(
    barcode = barcodes[ord],
    total.counts = as.numeric(total.counts[ord]),
    rank = seq_len(sum(keep)),
    stringsAsFactors = FALSE
  )
  
  if (cells > nrow(barcode.df)) {
    stop("cells exceeds the number of barcodes with positive counts.", call. = FALSE)
  }
  
  cutoff <- barcode.df$total.counts[cells]
  
  # ---- Create barcode rank plot ----
  graphics::plot(
    x = barcode.df$rank,
    y = barcode.df$total.counts,
    log = "xy",
    pch = 16,
    cex = 0.6,
    xlab = "Barcode rank",
    ylab = "Total counts",
    main = "Barcode-rank plot"
  )
  
  graphics::abline(v = cells, col = "red", lty = 2, lwd = 1.5)
  graphics::abline(h = cutoff, col = "red", lty = 2, lwd = 1.5)
  
  if (isTRUE(verbose)) {
    message("Selected barcode rank: ", cells)
    message("Count cutoff at this rank: ", cutoff)
  }
  
  attr(barcode.df, "selected.cells") <- cells
  attr(barcode.df, "cutoff") <- cutoff
  
  invisible(barcode.df)
}