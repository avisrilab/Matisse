# Matisse

<!-- badges: start -->
[![R-CMD-check](https://github.com/avisrilab/Matisse/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/avisrilab/Matisse/actions/workflows/R-CMD-check.yml)
[![pkgdown](https://github.com/avisrilab/Matisse/actions/workflows/pkgdown.yml/badge.svg)](https://avisrilab.github.io/Matisse)
<!-- badges: end -->

**Understand how genes are spliced differently across individual cells**

When a gene is expressed, its pre-mRNA can be spliced in multiple ways — different exons can be included or skipped, producing different protein isoforms. Matisse lets you measure this **alternative splicing** at single-cell resolution, so you can ask questions like:

- Do neurons splice *PTBP1* differently from astrocytes?
- Is the splicing change I see in bulk RNA-seq driven by one cell type or many?
- Which cells are switching isoforms along a developmental trajectory?

Matisse works alongside your existing [Seurat](https://satijalab.org/seurat/) workflow — your gene expression analysis stays intact, and splicing information is layered on top.

## What Matisse measures

For each cell and each splicing event, Matisse computes a **PSI value** (Percent Spliced In): the fraction of transcripts that *include* a given exon, ranging from 0 (exon always skipped) to 1 (exon always included). A PSI of 0.5 means half the transcripts in that cell include the exon.

## Installation

```r
# install.packages("remotes")
remotes::install_github("avisrilab/Matisse")
```

## Quick start

### Starting from junction counts (short-read 10x data)

If you ran STAR or STARsolo for alignment, you already have junction counts — a table recording how many reads span each exon-exon junction in each cell.

```r
library(Matisse)
library(Seurat)

# 1. Load your existing Seurat object (gene expression, UMAP, clusters, etc.)
seu <- readRDS("my_seurat.rds")

# 2. Load junction counts from STARsolo (cells × junctions)
jxn_counts <- readRDS("junction_counts.rds")

# 3. Load splice event definitions (from rMATS, SUPPA2, or BuildSimpleEvents)
event_data <- read.csv("events.csv")

# 4. Combine everything into a Matisse object
obj <- CreateMatisseObject(
  seurat          = seu,
  junction_counts = jxn_counts,
  event_data      = event_data
)

# 5. Calculate the PSI (splicing ratio) for each cell and each event
obj <- CalculatePSI(obj, min_coverage = 5)

# 6. Run quality control and remove low-quality cells and events
obj <- ComputeIsoformQC(obj)
obj <- FilterCells(obj, min_junctions = 5, min_pct_covered = 10)
obj <- FilterEvents(obj, min_cells_covered = 20)

# 7. Visualize — overlay splicing on your UMAP
PlotPSIUMAP(obj, event_id = "SE_PTBP1_e9")
PlotPSIViolin(obj, event_id = "SE_PTBP1_e9", group_by = "seurat_clusters")
```

### Starting from transcript counts (Bagpiper or long-read data)

If your quantifier reports counts per transcript isoform rather than per junction, use this path instead. Provide your transcript count table and SUPPA2 annotation files — Matisse handles the rest.

```r
# transcript_counts: transcripts × cells, with GENCODE transcript IDs as row names
obj <- CreateMatisseObjectFromTranscripts(
  seurat            = seu,
  transcript_counts = transcript_counts,
  ioe_files         = c("events_SE.ioe", "events_RI.ioe"),
  min_coverage      = 5L
)
```

From here, all QC, filtering, and visualization steps are identical to the junction-count workflow above.

## Documentation

Full vignettes and function reference are at **<https://avisrilab.github.io/Matisse>**.

## Citation

If you use Matisse in your research, please cite:

> Srivastava A. (2026). *Matisse: Multi-modal Analysis of Transcript Isoforms in Single-Cell Sequencing Experiments*. R package version 0.1.0. https://github.com/avisrilab/Matisse
