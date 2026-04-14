# Matisse 0.1.0

## Initial release

* New `MatisseObject` S4 class wrapping a `Seurat` object with isoform-resolved
  layers: junction counts, PSI matrix, inclusion/exclusion counts, event and
  junction annotation tables, and per-cell isoform metadata.
* `CreateMatisseObject()`: primary constructor with input validation and
  automatic cell-barcode alignment.
* `CalculatePSI()`: computes PSI matrices from junction count data and splice
  event definitions. Supports both `MatisseObject` and bare matrix input.
* `SummarizePSI()`: per-event summary statistics (mean, median, sd, coverage).
* `ComputeIsoformQC()`: per-cell QC metrics including junction detection rate,
  total junction reads, event coverage, and mean PSI.
* `FilterCells()` and `FilterEvents()`: threshold-based filtering with
  informative removal summaries.
* Visualization: `PlotPSIUMAP()`, `PlotPSIViolin()`, `PlotPSIHeatmap()`,
  `PlotJunctionCoverage()`, `PlotQCMetrics()`.
* `BuildSimpleEvents()`: convenience helper for generating one-vs-rest event
  tables from a junction annotation.
* `MergeMatisse()`: concatenate two `MatisseObject`s with cell-prefix deduplication.
* GitHub Actions workflows for R CMD check (multi-OS), pkgdown website
  deployment, and test coverage reporting.
