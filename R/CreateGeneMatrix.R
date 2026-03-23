#' Aggregate a transcript-level count matrix to a gene-level matrix
#'
#' @param tmat A sparse matrix with transcripts/isoforms as rows and cells as columns.
#' @param gtf.df A data.frame aligned to `tmat` with at least columns `transcript_id` and a gene label column (see ReadGTF()).
#' @param gene.col Character(1). Column name in `gtf.df` to use for grouping transcripts into genes (default = "gene").
#' @param verbose Logical(1). If TRUE, print progress messages (default = FALSE).
#'
#' @return A sparse gene-by-cell matrix (`Matrix::dgCMatrix`) with rownames = genes and colnames = barcodes.
#' @export
CreateGeneMatrix <- function(tmat,
                       gtf.df,
                       gene.col = "gene",
                       verbose = FALSE) {
  
  # ---- Verify parameter inputs ----
  if (!(inherits(tmat, "Matrix") || is.matrix(tmat))) {
    stop("tmat must be a matrix-like object with transcripts as rows and cells as columns.",
         call. = FALSE)
  }
  
  if (is.null(rownames(tmat)) || anyNA(rownames(tmat)) || any(!nzchar(rownames(tmat)))) {
    stop("tmat must have non-empty rownames containing transcript IDs.", call. = FALSE)
  }
  
  if (is.null(colnames(tmat)) || anyNA(colnames(tmat)) || any(!nzchar(colnames(tmat)))) {
    stop("tmat must have non-empty colnames containing cell barcodes.", call. = FALSE)
  }
  
  if (!inherits(tmat, "dgCMatrix")) {
    tmat <- methods::as(tmat, "dgCMatrix")
  }
  
  if (!is.data.frame(gtf.df)) {
    stop("gtf.df must be a data.frame.", call. = FALSE)
  }
  
  if (!("transcript_id" %in% names(gtf.df))) {
    stop("gtf.df must contain a 'transcript_id' column.", call. = FALSE)
  }
  
  if (!is.character(gene.col) || length(gene.col) != 1L || is.na(gene.col) || !nzchar(gene.col)) {
    stop("gene.col must be a non-empty character(1).", call. = FALSE)
  }
  
  if (!(gene.col %in% names(gtf.df))) {
    stop("gtf.df must contain the gene column: '", gene.col, "'.", call. = FALSE)
  }
  
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be a logical(1): TRUE or FALSE.", call. = FALSE)
  }
  
  if (nrow(gtf.df) != nrow(tmat)) {
    stop("gtf.df must have the same number of rows as tmat has transcripts (rows).",
         call. = FALSE)
  }
  
  tx.ids.mat <- trimws(as.character(rownames(tmat)))
  tx.ids.df  <- trimws(as.character(gtf.df$transcript_id))
  
  if (anyNA(tx.ids.df) || any(!nzchar(tx.ids.df))) {
    stop("gtf.df$transcript_id contains NA or empty values.", call. = FALSE)
  }
  
  if (!identical(tx.ids.df, tx.ids.mat)) {
    stop("gtf.df is not aligned to tmat: gtf.df$transcript_id must exactly match rownames(tmat) in the same order.",
         call. = FALSE)
  }
  
  if (anyDuplicated(tx.ids.mat)) {
    warning("Duplicate transcript IDs detected in tmat/gtf.df alignment.", call. = FALSE)
  }
  
  # ---- Extract gene labels ----
  gene <- trimws(as.character(gtf.df[[gene.col]]))
  missing.gene <- is.na(gene) | !nzchar(gene)
  n.missing <- sum(missing.gene)
  
  if (n.missing > 0L) {
    warning(n.missing,
            " transcripts have missing gene labels and will be excluded from aggregation.",
            call. = FALSE)
  }
  
  keep <- !missing.gene
  if (!any(keep)) {
    stop("No transcripts with non-missing gene labels are available for aggregation.",
         call. = FALSE)
  }
  
  tmat2 <- tmat[keep, , drop = FALSE]
  gene2 <- gene[keep]
  
  genes <- unique(gene2)
  
  if (isTRUE(verbose)) {
    message("Aggregating ", nrow(tmat2), " transcripts into ", length(genes), " genes.")
  }
  
  # ---- Create sparse transcript-to-gene mapping ----
  map <- Matrix::sparseMatrix(
    i = seq_along(gene2),
    j = match(gene2, genes),
    x = 1,
    dims = c(length(gene2), length(genes)),
    dimnames = list(rownames(tmat2), genes)
  )
  
  # t(map) %*% tmat2 gives gene x cell
  gmat <- Matrix::t(map) %*% tmat2
  
  if (!inherits(gmat, "dgCMatrix")) {
    gmat <- methods::as(gmat, "dgCMatrix")
  }
  
  rownames(gmat) <- genes
  colnames(gmat) <- colnames(tmat)
  
  gmat
}