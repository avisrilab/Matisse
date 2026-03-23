#' Create a Seurat object with transcript and gene assays plus transcript-level QC metrics
#'
#' @param tmat A sparse matrix with transcripts/isoforms as rows and cells as columns.
#' @param gmat A sparse matrix with genes as rows and cells as columns.
#' @param gtf.df A data.frame aligned to `tmat` with at least columns `transcript_id` and a gene label column.
#' @param project Character(1). Seurat project name (default = "SeuratProject").
#' @param transcript.assay Character(1). Name for transcript assay (default = "transcript")
#' @param gene.assay Character(1). Name for gene assay (default = "gene")
#' @param protein.coding.label Character(1). Label indicating protein-coding transcripts in `gtf.df$transcript.type` (default = "protein_coding").
#' @param mt.pattern Character(1). Regex used to identify mitochondrial genes in the gene assay. Default matches common human and mouse conventions (default = '^MT-|^mt-')"
#' @param verbose Logical(1). If TRUE, print progress messages (default = FALSE).
#'
#' @return A Seurat object with transcript and gene assays plus added QC metadata.
#' @export
CreateObject <- function(tmat,
                         gmat,
                         gtf.df,
                         project = "SeuratProject",
                         transcript.assay = "transcript",
                         gene.assay = "gene",
                         protein.coding.label = "protein_coding",
                         mt.pattern = "^MT-|^mt-",
                         verbose = FALSE) {
  
  # ---- Helper functions ----
  .check_scalar_string <- function(x, nm) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
      stop(nm, " must be a non-empty character(1).", call. = FALSE)
    }
  }
  
  .check_flag <- function(x, nm) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop(nm, " must be TRUE or FALSE.", call. = FALSE)
    }
  }
  
  # ---- Verify parameter inputs ----
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required. Please install it to create a Seurat object.", call. = FALSE)
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.", call. = FALSE)
  }
  
  if (!(inherits(tmat, "Matrix") || is.matrix(tmat))) {
    stop("tmat must be a matrix-like object with transcripts as rows and cells as columns.",
         call. = FALSE)
  }
  if (!(inherits(gmat, "Matrix") || is.matrix(gmat))) {
    stop("gmat must be a matrix-like object with genes as rows and cells as columns.",
         call. = FALSE)
  }
  
  if (is.null(rownames(tmat)) || anyNA(rownames(tmat)) || any(!nzchar(rownames(tmat)))) {
    stop("tmat must have non-empty rownames containing transcript IDs.", call. = FALSE)
  }
  if (is.null(colnames(tmat)) || anyNA(colnames(tmat)) || any(!nzchar(colnames(tmat)))) {
    stop("tmat must have non-empty colnames containing cell barcodes.", call. = FALSE)
  }
  
  if (is.null(rownames(gmat)) || anyNA(rownames(gmat)) || any(!nzchar(rownames(gmat)))) {
    stop("gmat must have non-empty rownames containing gene names.", call. = FALSE)
  }
  if (is.null(colnames(gmat)) || anyNA(colnames(gmat)) || any(!nzchar(colnames(gmat)))) {
    stop("gmat must have non-empty colnames containing cell barcodes.", call. = FALSE)
  }
  
  .check_scalar_string(project, "project")
  .check_scalar_string(transcript.assay, "transcript.assay")
  .check_scalar_string(gene.assay, "gene.assay")
  .check_scalar_string(protein.coding.label, "protein.coding.label")
  .check_scalar_string(mt.pattern, "mt.pattern")
  .check_flag(verbose, "verbose")
  
  if (identical(transcript.assay, gene.assay)) {
    stop("transcript.assay and gene.assay must be different.", call. = FALSE)
  }
  
  if (!identical(as.character(colnames(tmat)), as.character(colnames(gmat)))) {
    stop("tmat and gmat must have identical colnames (cells) in the same order.", call. = FALSE)
  }
  
  if (!is.data.frame(gtf.df)) {
    stop("gtf.df must be a data.frame.", call. = FALSE)
  }
  
  req.cols <- c("transcript_id", "gene", "transcript.type")
  missing.cols <- setdiff(req.cols, names(gtf.df))
  if (length(missing.cols) > 0L) {
    stop("gtf.df is missing required columns: ", paste(missing.cols, collapse = ", "),
         call. = FALSE)
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
  
  if (!inherits(tmat, "dgCMatrix")) tmat <- methods::as(tmat, "dgCMatrix")
  if (!inherits(gmat, "dgCMatrix")) gmat <- methods::as(gmat, "dgCMatrix")
  
  tx.type <- trimws(as.character(gtf.df$transcript.type))
  is.pc <- !is.na(tx.type) & nzchar(tx.type) &
    (tolower(tx.type) == tolower(protein.coding.label))
  
  if (sum(is.pc) == 0L) {
    warning("No transcripts matched protein.coding.label = '", protein.coding.label,
            "'. percent.protein.coding will be 0 for all cells.", call. = FALSE)
  }
  
  gene <- trimws(as.character(gtf.df$gene))
  missing.gene <- is.na(gene) | !nzchar(gene)
  
  if (sum(missing.gene) > 0L) {
    warning(sum(missing.gene),
            " transcripts have gene = NA/empty in gtf.df and will be excluded from avg.tx.per.gene and n.multiiso.genes.",
            call. = FALSE)
  }
  
  keep <- !missing.gene
  if (!any(keep)) {
    stop("No transcripts with non-missing gene labels available for transcript QC metrics.",
         call. = FALSE)
  }
  
  tmat.qc <- tmat[keep, , drop = FALSE]
  gene.qc <- gene[keep]
  
  f.gene <- factor(gene.qc)
  gene.int <- as.integer(f.gene)
  
  if (isTRUE(verbose)) {
    message("Computing transcript QC metrics using ",
            nrow(tmat.qc), " transcripts mapped to ",
            length(levels(f.gene)), " genes.")
  }
  
  # ---- Compute transcript-level QC metrics ----
  compute.basic.tx.qc <- function(tmat.csc, gene.int.vec) {
    nc <- ncol(tmat.csc)
    p <- tmat.csc@p
    i <- tmat.csc@i
    
    avg.tx.per.gene <- numeric(nc)
    n.multiiso.genes <- integer(nc)
    
    for (col in seq_len(nc)) {
      start <- p[col] + 1L
      end <- p[col + 1L]
      
      if (end < start) next
      
      rows <- i[start:end] + 1L
      gidx <- gene.int.vec[rows]
      
      n_tx_detected <- length(rows)
      ug <- unique(gidx)
      detected_genes <- length(ug)
      
      if (detected_genes == 0L) next
      
      avg.tx.per.gene[col] <- n_tx_detected / detected_genes
      
      gmap <- match(gidx, ug)
      iso_counts <- tabulate(gmap, nbins = length(ug))
      n.multiiso.genes[col] <- sum(iso_counts >= 2L)
    }
    
    list(
      avg.tx.per.gene = avg.tx.per.gene,
      n.multiiso.genes = n.multiiso.genes
    )
  }
  
  qc <- compute.basic.tx.qc(tmat.qc, gene.int)
  
  total.counts <- Matrix::colSums(tmat)
  pc.counts <- Matrix::colSums(tmat[is.pc, , drop = FALSE])
  percent.protein.coding <- ifelse(total.counts > 0, (pc.counts / total.counts) * 100, 0)
  
  n.type.missing <- sum(is.na(tx.type) | !nzchar(tx.type))
  if (n.type.missing > 0L) {
    warning(n.type.missing,
            " transcripts have missing transcript.type in gtf.df; they contribute to total counts but not protein-coding counts.",
            call. = FALSE)
  }
  
  if (isTRUE(verbose)) {
    message("Creating Seurat object with transcript assay '", transcript.assay, "'.")
  }
  
  # ---- Create Seurat object ----
  obj <- Seurat::CreateSeuratObject(
    counts = tmat,
    assay = transcript.assay,
    project = project
  )
  
  if (isTRUE(verbose)) {
    message("Adding gene assay '", gene.assay, "'.")
  }
  
  obj[[gene.assay]] <- Seurat::CreateAssayObject(counts = gmat)
  
  obj[["avg.tx.per.gene"]] <- qc$avg.tx.per.gene
  obj[["n.multiiso.genes"]] <- qc$n.multiiso.genes
  obj[["percent.protein.coding"]] <- as.numeric(percent.protein.coding)
  
  Seurat::DefaultAssay(obj) <- gene.assay
  obj <- Seurat::PercentageFeatureSet(obj, pattern = mt.pattern, col.name = "percent.mt")
  Seurat::DefaultAssay(obj) <- transcript.assay
  
  obj
}