# Matisse

**Multi-modal Analysis of Transcript Isoforms in Single-Cell Sequencing
Experiments**

Matisse is an R package for integrative analysis of isoform-resolved
single-cell transcriptomics. It extends
[Seurat](https://satijalab.org/seurat/) and
[Signac](https://stuartlab.org/signac/) with a dedicated `MatisseObject`
that co-stores gene expression and isoform layers in a single
synchronized object.

## Key capabilities

- **PSI matrix** — compute per-cell Percent Spliced In (PSI) values from
  junction read counts for any set of annotated splice events
- **Isoform QC** — per-cell metrics (junction detection rate, event
  coverage, mean PSI) with threshold-based filtering
- **Visualization** — UMAP overlays, violin plots, PSI heatmaps, and
  junction coverage plots via a consistent ggplot2 API
- **Seurat-native** — the embedded Seurat object is always in sync; all
  Seurat workflows continue to work unchanged

## Installation

``` r
# install.packages("remotes")
remotes::install_github("avisrilab/Matisse")
```

## Quick start

``` r
library(Matisse)
library(Seurat)

# 1. Start from an existing Seurat object
seu <- readRDS("my_seurat.rds")

# 2. Load junction counts (cells × junctions sparse matrix)
jxn_counts <- readRDS("junction_counts.rds")

# 3. Load splice event definitions (e.g. exported from rMATS / SUPPA2)
event_data <- read.csv("events.csv")

# 4. Build a MatisseObject
obj <- CreateMatisseObject(
  seurat        = seu,
  junction_counts = jxn_counts,
  event_data    = event_data
)

# 5. Calculate PSI
obj <- CalculatePSI(obj, min_coverage = 5)

# 6. QC and filter
obj <- ComputeIsoformQC(obj)
obj <- FilterCells(obj, min_junctions = 5, min_pct_covered = 10)
obj <- FilterEvents(obj, min_cells_covered = 20)

# 7. Visualize
PlotPSIUMAP(obj, event_id = "SE_PTBP1_e9")
PlotPSIViolin(obj, event_id = "SE_PTBP1_e9", group_by = "seurat_clusters")
```

## Documentation

Full documentation and vignettes are at
**<https://avisrilab.github.io/Matisse>**.

## Citation

If you use Matisse in your research, please cite:

> Srivastava A. (2026). *Matisse: Multi-modal Analysis of Transcript
> Isoforms in Single-Cell Sequencing Experiments*. R package version
> 0.1.0. <https://github.com/avisrilab/Matisse>
