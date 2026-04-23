# Script to generate pre-built vignette figures.
# Run once from the package root: Rscript vignettes/generate_figures.R
#
# Biology: mouse cortex single-cell RNA-seq, ~500 cells.
# Two populations — excitatory neurons and astrocytes — show a clear splicing
# switch at five well-characterised alternative exons.
# PTBP1 exon 9:  PSI ~0.08 in neurons (exon skipped), ~0.82 in astrocytes (included)
# NRXN1 SS4:     PSI ~0.72 in neurons (neuron-enriched inclusion), ~0.18 in astrocytes
# MAP2 exon 16:  PSI ~0.80 in neurons, ~0.28 in astrocytes
# FLNA exon 30:  PSI ~0.54 in both (no difference — negative control)
# MAPT exon 10:  PSI ~0.22 in neurons, ~0.18 in astrocytes (mild difference)

suppressPackageStartupMessages({
  pkgload::load_all(quiet = TRUE)   # load Matisse from source
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

set.seed(42)

# ---------------------------------------------------------------------------
# 1. Simulate cell populations
# ---------------------------------------------------------------------------
n_neurons <- 300
n_astro   <- 200
n_cells   <- n_neurons + n_astro
n_genes   <- 500

cells      <- paste0("Cell", seq_len(n_cells))
cell_types <- factor(c(rep("Neuron",    n_neurons),
                       rep("Astrocyte", n_astro)),
                     levels = c("Neuron", "Astrocyte"))
names(cell_types) <- cells

# Gene expression counts — add strong marker genes so Seurat clusters cleanly
gene_counts <- matrix(
  rpois(n_genes * n_cells, lambda = 5),
  nrow = n_genes, ncol = n_cells,
  dimnames = list(paste0("Gene", seq_len(n_genes)), cells)
)
# Neuron markers (Rbfox3/NeuN-like): Gene1–Gene5
gene_counts[1:5, seq_len(n_neurons)] <-
  matrix(rpois(5 * n_neurons, 60), nrow = 5)
# Astrocyte markers (Sox9/Aldh1l1-like): Gene6–Gene10
gene_counts[6:10, seq_len(n_astro) + n_neurons] <-
  matrix(rpois(5 * n_astro, 60), nrow = 5)

seu <- CreateSeuratObject(counts = gene_counts, min.cells = 0, min.features = 0)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, verbose = FALSE, npcs = 15)

# Manually set UMAP coordinates — two tight, well-separated clusters
umap_mat <- rbind(
  cbind(rnorm(n_neurons, mean = -3.0, sd = 0.7),
        rnorm(n_neurons, mean =  0.3, sd = 0.7)),   # neurons
  cbind(rnorm(n_astro,   mean =  3.0, sd = 0.7),
        rnorm(n_astro,   mean = -0.3, sd = 0.7))    # astrocytes
)
rownames(umap_mat) <- cells
colnames(umap_mat) <- c("UMAP_1", "UMAP_2")
seu[["umap"]]   <- SeuratObject::CreateDimReducObject(
  embeddings = umap_mat, key = "UMAP_", assay = "RNA")
seu$cell_type   <- cell_types
seu$seurat_clusters <- ifelse(cell_types == "Neuron", "0", "1")

# ---------------------------------------------------------------------------
# 2. Simulate junction counts for 5 SE splicing events (3 junctions each:
#    inc_up, inc_dn, exc).  Coordinates are simplified but realistic.
# ---------------------------------------------------------------------------
sim_event <- function(n, psi_target, depth_lambda = 35) {
  a   <- psi_target * 12
  b   <- (1 - psi_target) * 12
  psi <- rbeta(n, a, b)
  dep <- rpois(n, lambda = depth_lambda)
  inc <- as.integer(round(psi * dep))
  exc <- dep - inc
  list(inc = inc, exc = exc)
}

# Event definitions: SUPPA2 SE format — SE:chr:jxn1_start-jxn1_end:jxn2_start-jxn2_end:strand
# Each event has 3 junctions: inc_up (exon1→cassette), inc_dn (cassette→exon3), exc (skip)
events <- list(
  "SE:chr18:3433648-3434699:3434801-3436055:-" = list(gene = "PTBP1", chr = "chr18", strand = "-", neu = 0.08, ast = 0.82,
    jxn_up_s = 3433648L, jxn_up_e = 3434699L,
    jxn_dn_s = 3434801L, jxn_dn_e = 3436055L),
  "SE:chr17:90598735-90619999:90620100-90638554:+" = list(gene = "NRXN1", chr = "chr17", strand = "+", neu = 0.72, ast = 0.18,
    jxn_up_s = 90598735L, jxn_up_e = 90619999L,
    jxn_dn_s = 90620100L, jxn_dn_e = 90638554L),
  "SE:chr1:66049218-66057999:66058150-66070000:+"  = list(gene = "MAP2",  chr = "chr1",  strand = "+", neu = 0.80, ast = 0.28,
    jxn_up_s = 66049218L, jxn_up_e = 66057999L,
    jxn_dn_s = 66058150L, jxn_dn_e = 66070000L),
  "SE:chrX:153587380-153588499:153588650-153590000:+" = list(gene = "FLNA",  chr = "chrX",  strand = "+", neu = 0.54, ast = 0.52,
    jxn_up_s = 153587380L, jxn_up_e = 153588499L,
    jxn_dn_s = 153588650L, jxn_dn_e = 153590000L),
  "SE:chr11:77714194-77719999:77720150-77730000:-"  = list(gene = "MAPT",  chr = "chr11", strand = "-", neu = 0.22, ast = 0.18,
    jxn_up_s = 77714194L, jxn_up_e = 77719999L,
    jxn_dn_s = 77720150L, jxn_dn_e = 77730000L)
)
n_events <- length(events)
n_jxns   <- n_events * 3L  # inc_up, inc_dn, exc per event

jxn_ids <- as.vector(rbind(
  paste0(vapply(events, `[[`, character(1), "gene"), "_inc_up"),
  paste0(vapply(events, `[[`, character(1), "gene"), "_inc_dn"),
  paste0(vapply(events, `[[`, character(1), "gene"), "_exc")
))
jxn_mat <- matrix(0L, nrow = n_cells, ncol = n_jxns,
                  dimnames = list(cells, jxn_ids))

for (e in seq_along(events)) {
  ev      <- events[[e]]
  neu_sim <- sim_event(n_neurons, ev$neu)
  ast_sim <- sim_event(n_astro,   ev$ast)
  i_up  <- (e - 1L) * 3L + 1L
  i_dn  <- (e - 1L) * 3L + 2L
  i_exc <- (e - 1L) * 3L + 3L
  # Each inclusion junction carries ~half the inclusion reads
  jxn_mat[seq_len(n_neurons),           i_up]  <- as.integer(round(neu_sim$inc / 2))
  jxn_mat[seq_len(n_neurons),           i_dn]  <- as.integer(round(neu_sim$inc / 2))
  jxn_mat[seq_len(n_neurons),           i_exc] <- neu_sim$exc
  jxn_mat[seq_len(n_astro) + n_neurons, i_up]  <- as.integer(round(ast_sim$inc / 2))
  jxn_mat[seq_len(n_astro) + n_neurons, i_dn]  <- as.integer(round(ast_sim$inc / 2))
  jxn_mat[seq_len(n_astro) + n_neurons, i_exc] <- ast_sim$exc
}
jxn_sparse <- Matrix::Matrix(jxn_mat, sparse = TRUE)

# ---------------------------------------------------------------------------
# 3. Build event and junction annotation tables
# ---------------------------------------------------------------------------
event_df <- data.frame(
  event_id            = names(events),
  gene_id             = vapply(events, `[[`, character(1), "gene"),
  chr                 = vapply(events, `[[`, character(1), "chr"),
  strand              = vapply(events, `[[`, character(1), "strand"),
  event_type          = rep("SE", n_events),
  inclusion_junctions = paste0(
    vapply(events, `[[`, character(1), "gene"), "_inc_up", ";",
    vapply(events, `[[`, character(1), "gene"), "_inc_dn"),
  exclusion_junctions = paste0(
    vapply(events, `[[`, character(1), "gene"), "_exc"),
  stringsAsFactors    = FALSE
)

junction_df <- data.frame(
  junction_id = jxn_ids,
  chr         = rep(vapply(events, `[[`, character(1), "chr"),  each = 3L),
  start       = as.integer(c(rbind(
    vapply(events, `[[`, integer(1), "jxn_up_s"),
    vapply(events, `[[`, integer(1), "jxn_dn_s"),
    vapply(events, `[[`, integer(1), "jxn_up_s")   # exc shares upstream donor
  ))),
  end         = as.integer(c(rbind(
    vapply(events, `[[`, integer(1), "jxn_up_e"),
    vapply(events, `[[`, integer(1), "jxn_dn_e"),
    vapply(events, `[[`, integer(1), "jxn_dn_e")   # exc shares downstream acceptor
  ))),
  strand      = rep(vapply(events, `[[`, character(1), "strand"), each = 3L),
  gene_id     = rep(vapply(events, `[[`, character(1), "gene"),   each = 3L),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# 4. Create MatisseObject and run analysis
# ---------------------------------------------------------------------------
obj <- CreateMatisseObject(
  seurat          = seu,
  junction_counts = jxn_sparse,
  event_data      = event_df,
  junction_data   = junction_df,
  verbose         = FALSE
)
obj <- CalculatePSI(obj,       min_coverage = 5,  verbose = FALSE)
obj <- ComputeIsoformQC(obj,                      verbose = FALSE)
obj <- FilterCells(obj,        min_junctions = 2, verbose = FALSE)
obj <- FilterEvents(obj,       min_cells_covered = 10, verbose = FALSE)

# ---------------------------------------------------------------------------
# 5. Save figures
# ---------------------------------------------------------------------------
dir.create("vignettes/figures", showWarnings = FALSE, recursive = TRUE)

# Figure 1 — QC metrics by cell type
p1 <- PlotQCMetrics(obj, group_by = "cell_type")
ggsave("vignettes/figures/qc_metrics.png",   p1,
       width = 8, height = 5, dpi = 150, bg = "white")
message("Saved qc_metrics.png")

ptbp1_id <- "SE:chr18:3433648-3434699:3434801-3436055:-"
nrxn1_id <- "SE:chr17:90598735-90619999:90620100-90638554:+"

# Figure 2 — UMAP coloured by PTBP1 PSI
p2 <- PlotUMAP(obj,
  feature = ptbp1_id,
  pt_size = 1.2,
  title   = "PTBP1 exon 9 — PSI per cell")
ggsave("vignettes/figures/ptbp1_umap.png",   p2,
       width = 7, height = 5.5, dpi = 150, bg = "white")
message("Saved ptbp1_umap.png")

# Figure 3 — Violin plot comparing cell types
p3 <- PlotViolin(obj,
  feature  = ptbp1_id,
  group_by = "cell_type",
  title    = "PTBP1 exon 9 — PSI by cell type")
ggsave("vignettes/figures/ptbp1_violin.png", p3,
       width = 5, height = 5, dpi = 150, bg = "white")
message("Saved ptbp1_violin.png")

# Figure 4 — Heatmap of all five events ordered by cell type
p4 <- PlotHeatmap(obj, group_by = "cell_type", max_cells = 400)
ggsave("vignettes/figures/psi_heatmap.png",  p4,
       width = 7, height = 6, dpi = 150, bg = "white")
message("Saved psi_heatmap.png")

# Figure 5 — CoveragePlot: PTBP1 all cells pooled (junction mode)
p5 <- CoveragePlot(obj,
  event_id  = ptbp1_id,
  title     = "PTBP1 exon 9 — junction coverage")
ggsave("vignettes/figures/ptbp1_coverage.png", p5,
       width = 7, height = 4, dpi = 150, bg = "white")
message("Saved ptbp1_coverage.png")

# Figure 6 — CoveragePlot: PTBP1 faceted by cell type
p6 <- CoveragePlot(obj,
  event_id  = ptbp1_id,
  group_by  = "cell_type",
  title     = "PTBP1 exon 9 — coverage by cell type")
ggsave("vignettes/figures/ptbp1_coverage_by_type.png", p6,
       width = 10, height = 4, dpi = 150, bg = "white")
message("Saved ptbp1_coverage_by_type.png")

# Figure 7 — CoveragePlot: NRXN1 faceted by cell type
p7 <- CoveragePlot(obj,
  event_id  = nrxn1_id,
  group_by  = "cell_type",
  arc_scale = "sqrt",
  title     = "NRXN1 SS4 — coverage by cell type")
ggsave("vignettes/figures/nrxn1_coverage_by_type.png", p7,
       width = 10, height = 4, dpi = 150, bg = "white")
message("Saved nrxn1_coverage_by_type.png")

# Figure 8 — CoveragePlot: event mode (long-read simulation)
# Simulate a small long-read object for the event-mode CoveragePlot demo
lr_cells <- paste0("LRCell", seq_len(40L))
lr_types <- rep(c("TypeA", "TypeB"), each = 20L)
lr_genes <- paste0("Gene", seq_len(50L))
lr_expr  <- matrix(rpois(50L * 40L, 5L), nrow = 50L, ncol = 40L,
                   dimnames = list(lr_genes, lr_cells))
lr_seu   <- CreateSeuratObject(counts = lr_expr)
set.seed(99L)
lr_umap  <- matrix(c(rnorm(20L, -2, 0.6), rnorm(20L, 2, 0.6),
                     rnorm(20L,  0, 0.6), rnorm(20L, 0, 0.6)),
                   ncol = 2L, dimnames = list(lr_cells, c("UMAP_1", "UMAP_2")))
lr_seu[["umap"]]   <- SeuratObject::CreateDimReducObject(
  embeddings = lr_umap, key = "UMAP_")
lr_seu$cell_type   <- lr_types

lr_txs   <- c("tx_inc_a", "tx_inc_b", "tx_exc_a", "tx_exc_b")
lr_tx_mat <- matrix(0L, nrow = 4L, ncol = 40L,
                    dimnames = list(lr_txs, lr_cells))
set.seed(7L)
lr_tx_mat["tx_inc_a",  1:20]  <- sample(5L:12L, 20L, replace = TRUE)
lr_tx_mat["tx_inc_b",  1:20]  <- sample(4L:10L, 20L, replace = TRUE)
lr_tx_mat["tx_exc_a",  1:20]  <- sample(0L:2L,  20L, replace = TRUE)
lr_tx_mat["tx_exc_b",  1:20]  <- sample(0L:2L,  20L, replace = TRUE)
lr_tx_mat["tx_inc_a", 21:40]  <- sample(0L:2L,  20L, replace = TRUE)
lr_tx_mat["tx_inc_b", 21:40]  <- sample(0L:2L,  20L, replace = TRUE)
lr_tx_mat["tx_exc_a", 21:40]  <- sample(5L:12L, 20L, replace = TRUE)
lr_tx_mat["tx_exc_b", 21:40]  <- sample(4L:10L, 20L, replace = TRUE)

lr_ioe <- tempfile(fileext = ".ioe")
writeLines(c(
  "seqname\tgene_id\tinclusion_transcripts\ttotal_transcripts",
  paste(c("chr1", "ENSG1;SE:chr1:1201-2999:3201-4999:+",
          "tx_inc_a,tx_inc_b",
          "tx_inc_a,tx_inc_b,tx_exc_a,tx_exc_b"), collapse = "\t")
), lr_ioe)

lr_obj <- CreateMatisseObject(
  seurat            = lr_seu,
  transcript_counts = Matrix::Matrix(lr_tx_mat, sparse = TRUE),
  ioe_files         = lr_ioe,
  min_coverage      = 1L,
  verbose           = FALSE
)

p8 <- CoveragePlot(lr_obj,
  event_id  = "SE:chr1:1201-2999:3201-4999:+",
  group_by  = "cell_type",
  title     = "SE event — long-read coverage by cell type")
ggsave("vignettes/figures/lr_coverage_by_type.png", p8,
       width = 10, height = 4, dpi = 150, bg = "white")
message("Saved lr_coverage_by_type.png")

message("All figures written to vignettes/figures/")
