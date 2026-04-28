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

# Junction IDs encode coordinates so .parse_junction_names can derive
# chr/start/end/strand for sashimi plots without a separate junction_data
# argument (P7).
.junction_ids_fixture <- c(
  "chr1-1000-1500-+",
  "chr1-2000-2500-+",
  "chr1-3000-3500-+",
  "chr1-4000-4500-+",
  "chr1-5000-5500-+",
  "chr1-6000-6500-+"
)

make_junction_counts <- function(n_cells = 10L, n_jxns = 6L, seed = 42L) {
  set.seed(seed)
  cells <- paste0("Cell", seq_len(n_cells))
  jxns  <- .junction_ids_fixture[seq_len(n_jxns)]
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
    event_id           = c("SE-gene1-e2", "SE-gene1-e3"),
    gene_id            = c("gene1",        "gene1"),
    chr                = c("chr1",          "chr1"),
    strand             = c("+",             "+"),
    event_type         = c("SE",            "SE"),
    inclusion_features = c(paste(.junction_ids_fixture[1], .junction_ids_fixture[2], sep = ";"),
                            paste(.junction_ids_fixture[3], .junction_ids_fixture[4], sep = ";")),
    exclusion_features = c(.junction_ids_fixture[5], .junction_ids_fixture[6]),
    stringsAsFactors   = FALSE
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
  CreateMatisseObject(
    seurat          = make_seurat(),
    junction_counts = make_junction_counts(),
    events          = make_event_data(),
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
  CreateMatisseObject(
    seurat          = make_seurat_with_umap(),
    junction_counts = make_junction_counts(),
    events          = make_event_data(),
    min_coverage    = 1L,
    verbose         = FALSE
  )
}

# ---- MatisseObject from transcripts (transcript mode) ----------------------

make_matisse_from_transcripts <- function() {
  skip_if_not_installed("Seurat")
  CreateMatisseObject(
    seurat            = make_seurat(),
    transcript_counts = make_transcript_counts(),
    events            = make_ioe_file(),
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
# Junctions (encoded as chr-start-end-strand for auto-coord parsing):
#   jxn_up  : 1201 → 2999  (exon1 donor  → cassette acceptor)
#   jxn_dn  : 3201 → 4999  (cassette donor → exon3 acceptor)
#   jxn_exc : 1201 → 4999  (exon1 donor  → exon3 acceptor, skip)
#
# 20 cells: Cell1-10 = TypeA (high inclusion), Cell11-20 = TypeB (high exclusion)

.se_junction_ids <- c(
  "chr1-1201-2999-+",   # jxn_up   (event 1 inclusion)
  "chr1-3201-4999-+",   # jxn_dn   (event 1 inclusion)
  "chr1-1201-4999-+",   # jxn_exc  (events 1 & 2 exclusion)
  "chr1-5001-6999-+",   # jxn_up2  (event 2 inclusion)
  "chr1-7001-8999-+"    # jxn_dn2  (event 2 inclusion)
)

make_se_event_data <- function() {
  # Two SE events — required because Assay5 needs >=2 features (events).
  data.frame(
    event_id           = c("SE:chr1:1201-2999:3201-4999:+",
                           "SE:chr1:5001-6999:7001-8999:+"),
    gene_id            = rep("ENSG1", 2),
    chr                = rep("chr1", 2),
    strand             = rep("+", 2),
    event_type         = rep("SE", 2),
    inclusion_features = c(paste(.se_junction_ids[1], .se_junction_ids[2], sep = ";"),
                            paste(.se_junction_ids[4], .se_junction_ids[5], sep = ";")),
    exclusion_features = c(.se_junction_ids[3], .se_junction_ids[3]),
    stringsAsFactors   = FALSE
  )
}

make_se_junction_counts <- function(n_cells = 20L, seed = 42L) {
  set.seed(seed)
  cells <- paste0("Cell", seq_len(n_cells))
  typeA <- seq_len(n_cells / 2L)
  typeB <- seq(n_cells / 2L + 1L, n_cells)
  jxns  <- .se_junction_ids
  mat   <- matrix(0L, nrow = n_cells, ncol = length(jxns),
                  dimnames = list(cells, jxns))
  jxn_up   <- .se_junction_ids[1]
  jxn_dn   <- .se_junction_ids[2]
  jxn_exc  <- .se_junction_ids[3]
  jxn_up2  <- .se_junction_ids[4]
  jxn_dn2  <- .se_junction_ids[5]
  mat[typeA, jxn_up]  <- sample(8L:15L, length(typeA), replace = TRUE)
  mat[typeA, jxn_dn]  <- sample(8L:15L, length(typeA), replace = TRUE)
  mat[typeA, jxn_exc] <- sample(0L:3L,  length(typeA), replace = TRUE)
  mat[typeA, jxn_up2] <- sample(0L:3L,  length(typeA), replace = TRUE)
  mat[typeA, jxn_dn2] <- sample(0L:3L,  length(typeA), replace = TRUE)
  mat[typeB, jxn_up]  <- sample(0L:3L,  length(typeB), replace = TRUE)
  mat[typeB, jxn_dn]  <- sample(0L:3L,  length(typeB), replace = TRUE)
  mat[typeB, jxn_exc] <- sample(8L:15L, length(typeB), replace = TRUE)
  mat[typeB, jxn_up2] <- sample(8L:15L, length(typeB), replace = TRUE)
  mat[typeB, jxn_dn2] <- sample(8L:15L, length(typeB), replace = TRUE)
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
    events          = make_se_event_data(),
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
    events            = make_se_ioe_file(),
    min_coverage      = 1L,
    verbose           = FALSE
  )
}
