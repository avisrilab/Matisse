# Matisse: Multi-modal Analysis of Transcript Isoforms in Single-Cell Sequencing Experiments

Matisse provides an integrated framework for isoform-resolved
single-cell RNA-seq analysis, built on top of Seurat and Signac.

Key capabilities:

- **MatisseObject** – an S4 class that wraps a `Seurat` object and
  co-stores raw isoform counts, PSI matrices, and splice event
  annotations, keeping gene expression and isoform layers synchronised.

- **PSI calculation** –
  [`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
  computes per-cell Percent Spliced In values from raw junction or
  transcript counts. Works in both junction mode (STARsolo) and
  transcript mode (Bagpiper, FLAMES, LIQA).

- **Isoform QC** – per-cell metrics (`nCount_isoform`,
  `nFeature_isoform`, `nPercent_isoform`) are written automatically at
  construction and by
  [`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md);
  [`FilterCells`](https://avisrilab.github.io/Matisse/reference/FilterCells.md)
  and
  [`FilterEvents`](https://avisrilab.github.io/Matisse/reference/FilterEvents.md)
  enforce quality thresholds.

- **Visualization** – UMAP overlays, violin plots, PSI heatmaps, and
  sashimi junction-arc plots via a consistent ggplot2-based API.

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
