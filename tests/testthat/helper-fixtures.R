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
  suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
}

# ---- Full MatisseObject fixture (junction mode, no PSI yet) ----------------

make_matisse_object <- function() {
  # P5: CreateMatisseObject now folds PSI calc in by default. Most existing
  # tests assume the fixture returns a no-PSI object and call CalculatePSI
  # explicitly, so the fixture passes defer_psi = TRUE. Auto-PSI is
  # exercised by dedicated tests in test-MatisseObject.R.
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
    defer_psi       = TRUE,
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
  seu[["umap"]] <- suppressWarnings(SeuratObject::CreateDimReducObject(
    embeddings = coords, key = "UMAP_"
  ))
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
  CreateMatisseObject(
    seurat          = seu,
    junction_counts = jxn_mat,
    event_data      = ev_data,
    junction_data   = jd_data,
    min_coverage    = 1L,
    verbose         = FALSE
  )
}

# ---- MatisseObject from transcripts (transcript mode) ----------------------

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

# ---- Realistic SE event fixtures for PlotSashimi --------------------------
#
# Gene structure (chr1, + strand):
#   Exon 1:        900 – 1200
#   Cassette exon: 3000 – 3200
#   Exon 3:        5000 – 5300
#
# Junctions:
#   jxn_up  : 1201 → 2999  (exon1 donor  → cassette acceptor)
#   jxn_dn  : 3201 → 4999  (cassette donor → exon3 acceptor)
#   jxn_exc : 1201 → 4999  (exon1 donor  → exon3 acceptor, skip)
#
# 20 cells: Cell1-10 = TypeA (high inclusion), Cell11-20 = TypeB (high exclusion)

make_se_junction_data <- function() {
  # Two SE events sharing the upstream junction; second event uses jxn_up2/jxn_dn2.
  # Assay5 requires >=2 features, so we need >=2 events and >=2 junctions stored.
  data.frame(
    junction_id = c("jxn-up", "jxn-dn", "jxn-exc", "jxn-up2", "jxn-dn2"),
    chr         = rep("chr1", 5),
    start       = c(1201L, 3201L, 1201L, 5001L, 7001L),
    end         = c(2999L, 4999L, 4999L, 6999L, 8999L),
    strand      = rep("+", 5),
    gene_id     = rep("ENSG1", 5),
    stringsAsFactors = FALSE
  )
}

make_se_event_data <- function() {
  # Two SE events — required because Assay5 needs >=2 features (events).
  data.frame(
    event_id            = c("SE:chr1:1201-2999:3201-4999:+",
                            "SE:chr1:5001-6999:7001-8999:+"),
    gene_id             = rep("ENSG1", 2),
    chr                 = rep("chr1", 2),
    strand              = rep("+", 2),
    event_type          = rep("SE", 2),
    inclusion_junctions = c("jxn-up;jxn-dn", "jxn-up2;jxn-dn2"),
    exclusion_junctions = c("jxn-exc",        "jxn-exc"),
    stringsAsFactors    = FALSE
  )
}

make_se_junction_counts <- function(n_cells = 20L, seed = 42L) {
  set.seed(seed)
  cells <- paste0("Cell", seq_len(n_cells))
  typeA <- seq_len(n_cells / 2L)
  typeB <- seq(n_cells / 2L + 1L, n_cells)
  jxns  <- c("jxn-up", "jxn-dn", "jxn-exc", "jxn-up2", "jxn-dn2")
  mat   <- matrix(0L, nrow = n_cells, ncol = length(jxns),
                  dimnames = list(cells, jxns))
  mat[typeA, "jxn-up"]  <- sample(8L:15L, length(typeA), replace = TRUE)
  mat[typeA, "jxn-dn"]  <- sample(8L:15L, length(typeA), replace = TRUE)
  mat[typeA, "jxn-exc"] <- sample(0L:3L,  length(typeA), replace = TRUE)
  mat[typeA, "jxn-up2"] <- sample(0L:3L,  length(typeA), replace = TRUE)
  mat[typeA, "jxn-dn2"] <- sample(0L:3L,  length(typeA), replace = TRUE)
  mat[typeB, "jxn-up"]  <- sample(0L:3L,  length(typeB), replace = TRUE)
  mat[typeB, "jxn-dn"]  <- sample(0L:3L,  length(typeB), replace = TRUE)
  mat[typeB, "jxn-exc"] <- sample(8L:15L, length(typeB), replace = TRUE)
  mat[typeB, "jxn-up2"] <- sample(8L:15L, length(typeB), replace = TRUE)
  mat[typeB, "jxn-dn2"] <- sample(8L:15L, length(typeB), replace = TRUE)
  Matrix::Matrix(mat, sparse = TRUE)
}

# Short-read (junction mode): 20 cells, 2 types, SE event, PSI calculated
make_matisse_short_read <- function() {
  skip_if_not_installed("Seurat")
  n_cells <- 20L
  seu <- make_seurat(n_cells = n_cells)
  set.seed(7L)
  coords <- matrix(stats::rnorm(n_cells * 2L), nrow = n_cells, ncol = 2L,
                   dimnames = list(colnames(seu), c("UMAP_1", "UMAP_2")))
  seu[["umap"]]   <- suppressWarnings(SeuratObject::CreateDimReducObject(
    embeddings = coords, key = "UMAP_"))
  seu$cell_type   <- rep(c("TypeA", "TypeB"), each = n_cells / 2L)
  CreateMatisseObject(
    seurat          = seu,
    junction_counts = make_se_junction_counts(n_cells = n_cells),
    event_data      = make_se_event_data(),
    junction_data   = make_se_junction_data(),
    min_coverage    = 1L,
    verbose         = FALSE
  )
}

make_se_ioe_file <- function(file = tempfile(fileext = ".ioe")) {
  # Two events — Assay5 requires >=2 features (events).
  lines <- c(
    "seqname\tgene_id\tinclusion_transcripts\ttotal_transcripts",
    paste(c("chr1",
            "ENSG1;SE:chr1:1201-2999:3201-4999:+",
            "tx_inc_a,tx_inc_b",
            "tx_inc_a,tx_inc_b,tx_exc_a,tx_exc_b"),
          collapse = "\t"),
    paste(c("chr1",
            "ENSG1;SE:chr1:5001-6999:7001-8999:+",
            "tx_inc_c",
            "tx_inc_c,tx_exc_b"),
          collapse = "\t")
  )
  writeLines(lines, file)
  file
}

make_se_transcript_counts <- function(n_cells = 20L, seed = 42L) {
  set.seed(seed)
  cells <- paste0("Cell", seq_len(n_cells))
  typeA <- seq_len(n_cells / 2L)
  typeB <- seq(n_cells / 2L + 1L, n_cells)
  # tx_inc_c is the inclusion transcript for the second SE event
  txs   <- c("tx_inc_a", "tx_inc_b", "tx_exc_a", "tx_exc_b", "tx_inc_c")
  mat   <- matrix(0L, nrow = length(txs), ncol = n_cells,
                  dimnames = list(txs, cells))
  mat["tx_inc_a", typeA] <- sample(5L:12L, length(typeA), replace = TRUE)
  mat["tx_inc_b", typeA] <- sample(4L:10L, length(typeA), replace = TRUE)
  mat["tx_exc_a", typeA] <- sample(0L:2L,  length(typeA), replace = TRUE)
  mat["tx_exc_b", typeA] <- sample(0L:2L,  length(typeA), replace = TRUE)
  mat["tx_inc_c", typeA] <- sample(3L:8L,  length(typeA), replace = TRUE)
  mat["tx_inc_a", typeB] <- sample(0L:2L,  length(typeB), replace = TRUE)
  mat["tx_inc_b", typeB] <- sample(0L:2L,  length(typeB), replace = TRUE)
  mat["tx_exc_a", typeB] <- sample(5L:12L, length(typeB), replace = TRUE)
  mat["tx_exc_b", typeB] <- sample(4L:10L, length(typeB), replace = TRUE)
  mat["tx_inc_c", typeB] <- sample(0L:2L,  length(typeB), replace = TRUE)
  Matrix::Matrix(mat, sparse = TRUE)
}

# Long-read (transcript mode): 20 cells, 2 types, SE event via transcripts + IOE
make_matisse_long_read <- function() {
  skip_if_not_installed("Seurat")
  n_cells <- 20L
  seu <- make_seurat(n_cells = n_cells)
  set.seed(7L)
  coords <- matrix(stats::rnorm(n_cells * 2L), nrow = n_cells, ncol = 2L,
                   dimnames = list(colnames(seu), c("UMAP_1", "UMAP_2")))
  seu[["umap"]] <- suppressWarnings(SeuratObject::CreateDimReducObject(
    embeddings = coords, key = "UMAP_"))
  seu$cell_type <- rep(c("TypeA", "TypeB"), each = n_cells / 2L)
  CreateMatisseObject(
    seurat            = seu,
    transcript_counts = make_se_transcript_counts(n_cells = n_cells),
    ioe_files         = make_se_ioe_file(),
    min_coverage      = 1L,
    verbose           = FALSE
  )
}
