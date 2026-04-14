# Changelog

## Matisse 0.1.0

### Initial release

- New `MatisseObject` S4 class wrapping a `Seurat` object with
  isoform-resolved layers: junction counts, PSI matrix,
  inclusion/exclusion counts, event and junction annotation tables, and
  per-cell isoform metadata.
- [`CreateMatisseObject()`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md):
  primary constructor with input validation and automatic cell-barcode
  alignment.
- [`CalculatePSI()`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md):
  computes PSI matrices from junction count data and splice event
  definitions. Supports both `MatisseObject` and bare matrix input.
- [`SummarizePSI()`](https://avisrilab.github.io/Matisse/reference/SummarizePSI.md):
  per-event summary statistics (mean, median, sd, coverage).
- [`ComputeIsoformQC()`](https://avisrilab.github.io/Matisse/reference/ComputeIsoformQC.md):
  per-cell QC metrics including junction detection rate, total junction
  reads, event coverage, and mean PSI.
- [`FilterCells()`](https://avisrilab.github.io/Matisse/reference/FilterCells.md)
  and
  [`FilterEvents()`](https://avisrilab.github.io/Matisse/reference/FilterEvents.md):
  threshold-based filtering with informative removal summaries.
- Visualization:
  [`PlotPSIUMAP()`](https://avisrilab.github.io/Matisse/reference/PlotPSIUMAP.md),
  [`PlotPSIViolin()`](https://avisrilab.github.io/Matisse/reference/PlotPSIViolin.md),
  [`PlotPSIHeatmap()`](https://avisrilab.github.io/Matisse/reference/PlotPSIHeatmap.md),
  [`PlotJunctionCoverage()`](https://avisrilab.github.io/Matisse/reference/PlotJunctionCoverage.md),
  [`PlotQCMetrics()`](https://avisrilab.github.io/Matisse/reference/PlotQCMetrics.md).
- [`BuildSimpleEvents()`](https://avisrilab.github.io/Matisse/reference/BuildSimpleEvents.md):
  convenience helper for generating one-vs-rest event tables from a
  junction annotation.
- [`MergeMatisse()`](https://avisrilab.github.io/Matisse/reference/MergeMatisse.md):
  concatenate two `MatisseObject`s with cell-prefix deduplication.
- GitHub Actions workflows for R CMD check (multi-OS), pkgdown website
  deployment, and test coverage reporting.
