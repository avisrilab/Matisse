# Package index

## Creating a Matisse object

Start here. These functions combine your Seurat object with splicing
data into a single Matisse object ready for analysis. Use
[`CreateMatisseObject()`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md)
if you have junction counts from STAR/STARsolo, or
[`CreateMatisseObjectFromTranscripts()`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObjectFromTranscripts.md)
if you have isoform-level counts from Bagpiper or a long-read
quantifier.

- [`show(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
  [`dim(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
  [`` `[`( ``*`<MatisseObject>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
  [`` `[[`( ``*`<MatisseObject>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
  [`` `$`( ``*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
  : The MatisseObject S4 class
- [`CreateMatisseObject()`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md)
  : Create a MatisseObject
- [`CreateMatisseObjectFromTranscripts()`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObjectFromTranscripts.md)
  : Create a MatisseObject from transcript-level counts

## Retrieve or update your data

Functions for pulling specific data tables out of your Matisse object,
or putting updated tables back in.

- [`GetSeurat()`](https://avisrilab.github.io/Matisse/reference/GetSeurat.md)
  : Get the embedded Seurat object
- [`GetPSI()`](https://avisrilab.github.io/Matisse/reference/GetPSI.md)
  : Get the PSI matrix
- [`SetPSI()`](https://avisrilab.github.io/Matisse/reference/SetPSI.md)
  : Set the PSI matrix
- [`GetJunctionCounts()`](https://avisrilab.github.io/Matisse/reference/GetJunctionCounts.md)
  : Get raw junction count matrix
- [`GetInclusionCounts()`](https://avisrilab.github.io/Matisse/reference/GetInclusionCounts.md)
  : Get inclusion read count matrix
- [`GetExclusionCounts()`](https://avisrilab.github.io/Matisse/reference/GetExclusionCounts.md)
  : Get exclusion read count matrix
- [`GetEventData()`](https://avisrilab.github.io/Matisse/reference/GetEventData.md)
  : Get splice event annotation table
- [`GetJunctionData()`](https://avisrilab.github.io/Matisse/reference/GetJunctionData.md)
  : Get junction annotation table
- [`MatisseMeta()`](https://avisrilab.github.io/Matisse/reference/MatisseMeta.md)
  [`` `MatisseMeta<-`() ``](https://avisrilab.github.io/Matisse/reference/MatisseMeta.md)
  : Get or set cell-level isoform metadata
- [`AddIsoformMetadata()`](https://avisrilab.github.io/Matisse/reference/AddIsoformMetadata.md)
  : Add columns to the isoform metadata

## Calculate splicing ratios (PSI)

Compute PSI (Percent Spliced In) — the fraction of transcripts in each
cell that include a given exon. Values range from 0 (exon always
skipped) to 1 (exon always included).

- [`CalculatePSI()`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
  : Calculate PSI matrix from junction counts
- [`SummarizePSI()`](https://avisrilab.github.io/Matisse/reference/SummarizePSI.md)
  : Summarize PSI distribution across cells for each event

## Quality control and filtering

Identify and remove cells or splicing events that don’t have enough data
for reliable analysis.

- [`ComputeIsoformQC()`](https://avisrilab.github.io/Matisse/reference/ComputeIsoformQC.md)
  : Compute per-cell isoform QC metrics
- [`FilterCells()`](https://avisrilab.github.io/Matisse/reference/FilterCells.md)
  : Filter cells by isoform QC thresholds
- [`FilterEvents()`](https://avisrilab.github.io/Matisse/reference/FilterEvents.md)
  : Filter splice events by coverage or variance

## Visualisation

Plot splicing patterns across your cells. Overlay PSI values on a UMAP,
compare splicing between cell types, or inspect junction usage for a
gene of interest.

- [`PlotPSIUMAP()`](https://avisrilab.github.io/Matisse/reference/PlotPSIUMAP.md)
  : UMAP plot colored by PSI of a specific splice event
- [`PlotPSIViolin()`](https://avisrilab.github.io/Matisse/reference/PlotPSIViolin.md)
  : Violin plot of PSI values split by cell group
- [`PlotPSIHeatmap()`](https://avisrilab.github.io/Matisse/reference/PlotPSIHeatmap.md)
  : Heatmap of PSI values across cells and events
- [`PlotJunctionCoverage()`](https://avisrilab.github.io/Matisse/reference/PlotJunctionCoverage.md)
  : Per-cell junction coverage plot for a gene
- [`PlotQCMetrics()`](https://avisrilab.github.io/Matisse/reference/PlotQCMetrics.md)
  : Violin/ridge plot of isoform QC metrics

## Utilities

Helper functions for building event tables and combining datasets.

- [`BuildSimpleEvents()`](https://avisrilab.github.io/Matisse/reference/BuildSimpleEvents.md)
  : Build a minimal junction event annotation table
- [`MergeMatisse()`](https://avisrilab.github.io/Matisse/reference/MergeMatisse.md)
  : Merge two MatisseObjects by cells
