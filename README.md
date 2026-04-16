# Matisse

<!-- badges: start -->
[![R-CMD-check](https://github.com/avisrilab/Matisse/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/avisrilab/Matisse/actions/workflows/R-CMD-check.yml)
[![pkgdown](https://github.com/avisrilab/Matisse/actions/workflows/pkgdown.yml/badge.svg)](https://avisrilab.github.io/Matisse)
<!-- badges: end -->

**Understand your cells, layer by layer — splicing and chromatin together**

Most single-cell tools stop at gene expression. Matisse goes further: it measures *which isoform* each cell is making and links those choices to chromatin accessibility in the same cells. Some cell types skip exons. Others include them. Open chromatin at a splice site may predict which. Matisse lets you see all of this, in the same object, on the same UMAP.

![PSI values for PTBP1 exon 9 on a UMAP of mouse cortex cells. Neurons (left) consistently skip this exon; astrocytes (right) include it.](man/figures/ptbp1_umap.png)

*Neurons and astrocytes express the same gene — but splice it differently. Matisse finds these differences automatically, alongside your existing chromatin and gene expression layers.*

---

## What you can discover

- **Cell-type-specific splicing** — Do my neurons and astrocytes process this exon differently?
- **Chromatin shapes isoforms** — Is the splicing switch I see linked to chromatin accessibility changes at the same locus?
- **Splicing along a trajectory** — Is there a coordinated isoform switch as cells differentiate?
- **Context for bulk data** — I see a splicing change in bulk RNA-seq — which cell type is driving it?

---

## How it works

For each cell and each splicing event, Matisse calculates a **PSI value** (Percent Spliced In): the fraction of that cell's transcripts that include a given exon.

- **PSI = 1** — every transcript in that cell includes the exon
- **PSI = 0** — every transcript skips the exon
- **PSI = 0.5** — half include, half skip

PSI values are stored alongside your gene expression and chromatin data in a single object, so you can overlay splicing on any plot you've already made — UMAPs, violin plots, heatmaps — and directly compare with ATAC peaks in the same cells.

---

## Works with your existing setup

Matisse layers on top of [Seurat](https://satijalab.org/seurat/) and [Signac](https://stuartlab.org/signac/) — your clustering, UMAP, chromatin peaks, and cell-type labels stay exactly as they are.

**Short-read RNA data (10x Chromium):**
Aligned with STAR or STARsolo → junction count matrix → Matisse

**Long-read / isoform-resolved data:**
Quantified with Bagpiper, FLAMES, or LIQA → transcript count table + SUPPA2 `.ioe` files → Matisse

**ATAC / chromatin accessibility (10x Multiome, Signac, ArchR):**
Your existing Signac object is embedded directly — chromatin and splicing layers stay synchronized

**Splice event annotations:**
Compatible with SUPPA2 `generateEvents`, rMATS output, or a simple table you build yourself.

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("avisrilab/Matisse")
```

---

## Quick start

### Starting from junction counts (10x / STARsolo)

```r
library(Matisse)
library(Seurat)

# 1. Your existing Seurat object — clusters, UMAP, everything intact
seu <- readRDS("my_seurat.rds")

# 2. Junction count matrix from STARsolo (cells × junctions)
jxn_counts <- readRDS("junction_counts.rds")

# 3. Splice event table (from SUPPA2, rMATS, or BuildSimpleEvents)
event_data <- read.csv("events.csv")

# 4. Build the Matisse object
obj <- CreateMatisseObject(
  seurat          = seu,
  junction_counts = jxn_counts,
  event_data      = event_data
)

# 5. Calculate PSI for every cell and every event
obj <- CalculatePSI(obj, min_coverage = 5)

# 6. Quality control
obj <- ComputeIsoformQC(obj)
obj <- FilterCells(obj, min_junctions = 5, min_pct_covered = 10)
obj <- FilterEvents(obj, min_cells_covered = 20)

# 7. Visualise — overlay splicing on your UMAP
PlotPSIUMAP(obj, event_id = "SE_PTBP1_e9")
PlotPSIViolin(obj, event_id = "SE_PTBP1_e9", group_by = "seurat_clusters")
```

### Starting from transcript counts (Bagpiper / long-read)

```r
obj <- CreateMatisseObjectFromTranscripts(
  seurat            = seu,
  transcript_counts = transcript_counts,
  ioe_files         = c("events_SE.ioe", "events_RI.ioe"),
  min_coverage      = 5L
)
```

From here, QC, filtering, and visualisation are identical to the junction-count workflow.

---

## Documentation

Full walkthrough and function reference: **<https://avisrilab.github.io/Matisse>**

---

## Citation

If you use Matisse in your research, please cite:

> Srivastava A. (2026). *Matisse: Multi-modal Analysis of Transcript Isoforms in Single-Cell Sequencing Experiments*. R package version 0.1.0. https://github.com/avisrilab/Matisse
