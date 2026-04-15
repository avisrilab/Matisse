# ---------------------------------------------------------------------------
# Shared test fixtures
# ---------------------------------------------------------------------------

# ---- Transcript count matrix (transcripts x cells) -------------------------

#' Create a tiny transcript count matrix.
#' Transcripts: tx1–tx8. Cells: Cell1–Cell10.
#' tx1+tx2 = inclusion for event 1, tx3 = exclusion for event 1
#' tx4+tx5 = inclusion for event 2, tx6+tx7 = exclusion for event 2
#' tx8     = not part of any event (should be ignored)
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

#' Write a SUPPA2 IOE file to a temp path and return the path.
#' Two events:
#'   SE_gene1_t12: inclusion = tx1,tx2  exclusion = tx3
#'   SE_gene1_t45: inclusion = tx4,tx5  exclusion = tx6,tx7
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
# These helpers are auto-loaded by testthat before test files run.
# They create minimal, self-contained objects without requiring real data.

# ---- Junction count matrix -------------------------------------------------

#' Create a tiny junction count matrix (cells x junctions)
make_junction_counts <- function(n_cells = 10L, n_jxns = 6L, seed = 42L) {
  set.seed(seed)
  cells <- paste0("Cell", seq_len(n_cells))
  jxns  <- paste0("jxn", seq_len(n_jxns))

  # Simulate sparse-ish counts
  mat <- matrix(
    sample(c(0L, 0L, 0L, 1L:20L), n_cells * n_jxns, replace = TRUE),
    nrow = n_cells, ncol = n_jxns,
    dimnames = list(cells, jxns)
  )
  Matrix::Matrix(mat, sparse = TRUE)
}

# ---- Event annotation table ------------------------------------------------

#' Two simple SE-type events using jxn1–jxn4
make_event_data <- function() {
  data.frame(
    event_id             = c("SE_gene1_e2", "SE_gene1_e3"),
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

# ---- Seurat object (created only when Seurat is installed) -----------------

#' Build a minimal Seurat object with 10 cells and 20 genes.
#' Skips the calling test if Seurat is not installed.
make_seurat <- function(n_cells = 10L, n_genes = 20L, seed = 1L) {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")

  set.seed(seed)
  cells <- paste0("Cell", seq_len(n_cells))
  genes <- paste0("Gene", seq_len(n_genes))

  counts <- matrix(
    rpois(n_genes * n_cells, lambda = 5),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(genes, cells)
  )

  Seurat::CreateSeuratObject(counts = counts)
}

# ---- Full MatisseObject fixture --------------------------------------------

#' Build a complete MatisseObject fixture with junctions, events, and QC.
make_matisse_object <- function() {
  skip_if_not_installed("Seurat")
  seu     <- make_seurat()
  jxn_mat <- make_junction_counts()
  ev_data <- make_event_data()
  jd_data <- make_junction_data()

  CreateMatisseObject(
    seurat        = seu,
    junction_counts = jxn_mat,
    event_data    = ev_data,
    junction_data = jd_data,
    verbose       = FALSE
  )
}
