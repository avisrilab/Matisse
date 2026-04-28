# Matisse: Multi-modal Analysis of Transcript Isoforms in Single-Cell Sequencing Experiments

Matisse provides an integrated framework for isoform-resolved
single-cell RNA-seq analysis, built on top of Seurat and Signac.

Key capabilities:

- **MatisseObject** – an S4 class that wraps a `Seurat` object and
  co-stores raw isoform counts, PSI matrices, and splice event
  annotations, keeping gene expression and isoform layers synchronised.

- **Single-step construction with PSI** –
  [`CreateMatisseObject`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md)
  computes per-cell Percent Spliced In values as part of construction
  (pass `defer_psi = TRUE` to skip). Works in both junction mode
  (STARsolo) and transcript mode (Bagpiper, FLAMES, LIQA).

- **Isoform QC** – `nCount_isoform` and `nFeature_isoform` are written
  at construction; `nPercent_isoform` is written by the PSI step.
  [`FilterCells`](https://avisrilab.github.io/Matisse/reference/FilterCells.md)
  and
  [`FilterEvents`](https://avisrilab.github.io/Matisse/reference/FilterEvents.md)
  enforce quality thresholds.

- **Visualization** – UMAP overlays, violin plots, PSI heatmaps, and
  sashimi junction-arc plots via a consistent ggplot2-based API.

## Storage convention

Matisse stores everything inside the embedded Seurat object. To access
raw data with native Seurat verbs:

- Junction or transcript counts:

  `SeuratObject::GetAssayData(GetSeurat(obj)[["isoform"]], "counts")`

- PSI values (events x cells):

  `SeuratObject::GetAssayData(GetSeurat(obj)[["psi"]], "data")`

- Inclusion / exclusion read counts:

  `SeuratObject::GetAssayData(GetSeurat(obj)[["psi"]], "counts")` and
  `... ([["psi"]], "exclusion")`

- Per-event annotation:

  `GetSeurat(obj)[["psi"]][[]]` (data.frame keyed by event_id rownames)

- Per-junction coordinates (junction mode):

  `GetSeurat(obj)[["isoform"]][[]]` – chr/start/end/strand auto-parsed
  from junction IDs at construction.

- Per-cell metadata:

  `MatisseMeta(obj)` or `GetSeurat(obj)@meta.data`

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
