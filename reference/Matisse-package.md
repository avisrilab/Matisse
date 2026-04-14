# Matisse: Multi-modal Analysis of Transcript Isoforms in Single-Cell Sequencing Experiments

Matisse provides an integrated framework for isoform-resolved
single-cell RNA-seq analysis, built on top of Seurat and Signac.

Key capabilities:

- **MatisseObject** — an S4 class that wraps a `Seurat` object and
  co-stores junction counts, PSI matrices, and splice event annotations,
  keeping gene expression and isoform layers synchronised.

- **PSI calculation** —
  [`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
  computes per-cell Percent Spliced In values from raw junction read
  counts and a user-supplied or auto-generated event annotation table.

- **Isoform QC** —
  [`ComputeIsoformQC`](https://avisrilab.github.io/Matisse/reference/ComputeIsoformQC.md)
  derives per-cell metrics (junction detection rate, event coverage,
  mean PSI);
  [`FilterCells`](https://avisrilab.github.io/Matisse/reference/FilterCells.md)
  and
  [`FilterEvents`](https://avisrilab.github.io/Matisse/reference/FilterEvents.md)
  enforce quality thresholds.

- **Visualization** — UMAP overlays, violin plots, PSI heatmaps, and
  junction coverage bar charts via a consistent ggplot2-based API.

## Package website

Full documentation and vignettes are available at
<https://avisrilab.github.io/Matisse>.

## See also

Useful links:

- <https://avisrilab.github.io/Matisse>

- <https://github.com/avisrilab/Matisse>

- Report bugs at <https://github.com/avisrilab/Matisse/issues>

## Author

**Maintainer**: k3yavi <avisrilab@gmail.com>
