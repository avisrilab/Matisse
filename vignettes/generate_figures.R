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
# 2. Simulate junction counts for 5 splicing events (2 junctions each)
# ---------------------------------------------------------------------------
sim_event <- function(n, psi_target, depth_lambda = 35) {
  # Draw PSI from a beta distribution centred on psi_target
  a   <- psi_target * 12
  b   <- (1 - psi_target) * 12
  psi <- rbeta(n, a, b)
  dep <- rpois(n, lambda = depth_lambda)
  inc <- as.integer(round(psi * dep))
  exc <- dep - inc
  list(inc = inc, exc = exc)
}

events <- list(
  "PTBP1:SE:chr18:3433647-3436055"     = list(neu = 0.08, ast = 0.82),
  "NRXN1:SS4:chr17:90598735-90638554"  = list(neu = 0.72, ast = 0.18),
  "MAP2:SE:chr1:66049218-66070000"     = list(neu = 0.80, ast = 0.28),
  "FLNA:SE:chrX:153587380-153590000"   = list(neu = 0.54, ast = 0.52),
  "MAPT:SE:chr11:77714194-77730000"    = list(neu = 0.22, ast = 0.18)
)
n_events <- length(events)
n_jxns   <- n_events * 2

jxn_mat <- matrix(0L, nrow = n_cells, ncol = n_jxns,
                  dimnames = list(cells, paste0("jxn", seq_len(n_jxns))))

for (e in seq_along(events)) {
  ev      <- events[[e]]
  neu_sim <- sim_event(n_neurons, ev$neu)
  ast_sim <- sim_event(n_astro,   ev$ast)
  inc_col <- (e - 1L) * 2L + 1L
  exc_col <- (e - 1L) * 2L + 2L
  jxn_mat[seq_len(n_neurons),          inc_col] <- neu_sim$inc
  jxn_mat[seq_len(n_neurons),          exc_col] <- neu_sim$exc
  jxn_mat[seq_len(n_astro) + n_neurons, inc_col] <- ast_sim$inc
  jxn_mat[seq_len(n_astro) + n_neurons, exc_col] <- ast_sim$exc
}
jxn_sparse <- Matrix::Matrix(jxn_mat, sparse = TRUE)

# ---------------------------------------------------------------------------
# 3. Build event annotation table
# ---------------------------------------------------------------------------
event_df <- data.frame(
  event_id            = names(events),
  gene_id             = c("PTBP1", "NRXN1", "MAP2", "FLNA", "MAPT"),
  chr                 = c("chr18", "chr17", "chr1", "chrX", "chr11"),
  strand              = c("-",     "+",     "+",    "+",    "-"),
  event_type          = rep("SE", n_events),
  inclusion_junctions = paste0("jxn", seq(1L, n_jxns - 1L, by = 2L)),
  exclusion_junctions = paste0("jxn", seq(2L, n_jxns,      by = 2L)),
  stringsAsFactors    = FALSE
)

# ---------------------------------------------------------------------------
# 4. Create MatisseObject and run analysis
# ---------------------------------------------------------------------------
obj <- CreateMatisseObject(
  seurat          = seu,
  junction_counts = jxn_sparse,
  event_data      = event_df,
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

# Figure 2 — UMAP coloured by PTBP1 PSI
p2 <- PlotPSIUMAP(obj,
  event_id  = "PTBP1:SE:chr18:3433647-3436055",
  pt_size   = 1.2,
  title     = "PTBP1 exon 9 — PSI per cell")
ggsave("vignettes/figures/ptbp1_umap.png",   p2,
       width = 7, height = 5.5, dpi = 150, bg = "white")
message("Saved ptbp1_umap.png")

# Figure 3 — Violin plot comparing cell types
p3 <- PlotPSIViolin(obj,
  event_id  = "PTBP1:SE:chr18:3433647-3436055",
  group_by  = "cell_type",
  title     = "PTBP1 exon 9 — PSI by cell type")
ggsave("vignettes/figures/ptbp1_violin.png", p3,
       width = 5, height = 5, dpi = 150, bg = "white")
message("Saved ptbp1_violin.png")

# Figure 4 — Heatmap of all five events ordered by cell type
p4 <- PlotPSIHeatmap(obj, group_by = "cell_type", max_cells = 400)
ggsave("vignettes/figures/psi_heatmap.png",  p4,
       width = 7, height = 6, dpi = 150, bg = "white")
message("Saved psi_heatmap.png")

message("All figures written to vignettes/figures/")
