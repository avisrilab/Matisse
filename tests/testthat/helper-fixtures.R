# ---------------------------------------------------------------------------
# Shared test fixtures
# Auto-loaded by testthat before test files run.
# ---------------------------------------------------------------------------

# ---- Transcript count matrix (transcripts x cells) -------------------------

make_transcript_counts <- function(n_cells = 10L, seed = 42L) {
  set.seed(seed)
  cells <- paste0("Cell", seq_len(n_cells))
  txs   <- paste0("tx",   seq_len(8L))
  mat   <- matrix(
    sample(c(0L, 0L, 1L:15L), n_cells * 8L, replace = TRUE),
    nrow = 8L, ncol = n_cells,
    dimnames = list(txs, cells)
  )
  Matrix::Matrix(mat, sparse = TRUE)
}

#' Write a SUPPA2 IOE file (two SE events).
make_ioe_file <- function(file = tempfile(fileext = ".ioe")) {
  lines <- c(
    "seqname\tgene_id\tinclusion_transcripts\ttotal_transcripts",
    paste(c("chr1",
            "ENSG00000001;SE:chr1:100-200:300-400:+",
            "tx1,tx2",
            "tx1,tx2,tx3"),
          collapse = "\t"),
    paste(c("chr1",
            "ENSG00000001;SE:chr1:500-600:700-800:+",
            "tx4,tx5",
            "tx4,tx5,tx6,tx7"),
          collapse = "\t")
  )
  writeLines(lines, file)
  file
}

# ---- Junction count matrix (cells x junctions) -----------------------------

make_junction_counts <- function(n_cells = 10L, n_jxns = 6L, seed = 42L) {
  set.seed(seed)
  cells <- paste0("Cell", seq_len(n_cells))
  jxns  <- paste0("jxn",  seq_len(n_jxns))
  mat   <- matrix(
    sample(c(0L, 0L, 0L, 1L:20L), n_cells * n_jxns, replace = TRUE),
    nrow = n_cells, ncol = n_jxns,
    dimnames = list(cells, jxns)
  )
  Matrix::Matrix(mat, sparse = TRUE)
}

# ---- Event annotation table ------------------------------------------------

make_event_data <- function() {
  data.frame(
    event_id             = c("SE-gene1-e2", "SE-gene1-e3"),
    gene_id              = c("gene1",        "gene1"),
    chr                  = c("chr1",          "chr1"),
    strand               = c("+",             "+"),
    event_type           = c("SE",            "SE"),
    inclusion_junctions  = c("jxn1;jxn2",    "jxn3;jxn4"),
    exclusion_junctions  = c("jxn5",          "jxn6"),
    stringsAsFactors     = FALSE
  )
}

# ---- Junction annotation table ---------------------------------------------

make_junction_data <- function() {
  data.frame(
    junction_id = paste0("jxn", 1:6),
    chr         = rep("chr1", 6),
    start       = c(1000L, 2000L, 3000L, 4000L, 5000L, 6000L),
    end         = c(1500L, 2500L, 3500L, 4500L, 5500L, 6500L),
    strand      = rep("+", 6),
    gene_id     = rep("gene1", 6),
    stringsAsFactors = FALSE
  )
}

# ---- Seurat object ---------------------------------------------------------

make_seurat <- function(n_cells = 10L, n_genes = 20L, seed = 1L) {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  set.seed(seed)
  cells  <- paste0("Cell", seq_len(n_cells))
  genes  <- paste0("Gene", seq_len(n_genes))
  counts <- matrix(rpois(n_genes * n_cells, lambda = 5),
                   nrow = n_genes, ncol = n_cells,
                   dimnames = list(genes, cells))
  Seurat::CreateSeuratObject(counts = counts)
}

# ---- Full MatisseObject fixture (junction mode, no PSI yet) ----------------

make_matisse_object <- function() {
  skip_if_not_installed("Seurat")
  seu     <- make_seurat()
  jxn_mat <- make_junction_counts()
  ev_data <- make_event_data()
  jd_data <- make_junction_data()
  CreateMatisseObject(
    seurat          = seu,
    junction_counts = jxn_mat,
    event_data      = ev_data,
    junction_data   = jd_data,
    verbose         = FALSE
  )
}

# ---- Seurat object with UMAP and cell_type metadata -----------------------

make_seurat_with_umap <- function(n_cells = 10L, n_genes = 20L, seed = 1L) {
  seu <- make_seurat(n_cells = n_cells, n_genes = n_genes, seed = seed)
  set.seed(seed + 100L)
  coords <- matrix(
    stats::rnorm(n_cells * 2L),
    nrow     = n_cells,
    ncol     = 2L,
    dimnames = list(colnames(seu), c("UMAP_1", "UMAP_2"))
  )
  seu[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = coords, key = "UMAP_"
  )
  seu$cell_type       <- rep(c("TypeA", "TypeB"), length.out = n_cells)
  seu$seurat_clusters <- factor(rep(0L, n_cells))
  seu
}

# ---- MatisseObject with UMAP and PSI calculated ----------------------------

make_matisse_with_umap <- function() {
  skip_if_not_installed("Seurat")
  seu     <- make_seurat_with_umap()
  jxn_mat <- make_junction_counts()
  ev_data <- make_event_data()
  jd_data <- make_junction_data()
  obj <- CreateMatisseObject(
    seurat          = seu,
    junction_counts = jxn_mat,
    event_data      = ev_data,
    junction_data   = jd_data,
    verbose         = FALSE
  )
  CalculatePSI(obj, min_coverage = 1L, verbose = FALSE)
}

# ---- MatisseObject with QC computed ----------------------------------------

make_matisse_with_qc <- function() {
  obj <- make_matisse_with_umap()
  ComputeIsoformQC(obj, verbose = FALSE)
}

# ---- MatisseObject from transcripts (event mode) ---------------------------

make_matisse_from_transcripts <- function() {
  skip_if_not_installed("Seurat")
  seu    <- make_seurat()
  tx_mat <- make_transcript_counts()
  f      <- make_ioe_file()
  CreateMatisseObject(
    seurat            = seu,
    transcript_counts = tx_mat,
    ioe_files         = f,
    min_coverage      = 1L,
    verbose           = FALSE
  )
}
