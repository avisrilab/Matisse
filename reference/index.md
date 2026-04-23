# Package index

## Creating a Matisse object

Start here. This function combines your Seurat object with splicing data
into a single Matisse object ready for analysis. Pass `junction_counts`
for short-read (junction) mode, or `transcript_counts` + `ioe_files` for
long-read (event) mode.

- [`show(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
  [`dim(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
  [`` `[`( ``*`<MatisseObject>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
  [`` `[[`( ``*`<MatisseObject>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
  [`` `$`( ``*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
  : The MatisseObject S4 class
- [`CreateMatisseObject()`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md)
  : Create a MatisseObject

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
- [`GetTranscriptCounts()`](https://avisrilab.github.io/Matisse/reference/GetTranscriptCounts.md)
  : Get transcript count matrix
- [`GetEventData()`](https://avisrilab.github.io/Matisse/reference/GetEventData.md)
  : Get splice event annotation table
- [`GetJunctionData()`](https://avisrilab.github.io/Matisse/reference/GetJunctionData.md)
  : Get junction annotation table
- [`MatisseMeta()`](https://avisrilab.github.io/Matisse/reference/MatisseMeta.md)
  [`` `MatisseMeta<-`() ``](https://avisrilab.github.io/Matisse/reference/MatisseMeta.md)
  : Get or set cell-level metadata
- [`AddIsoformMetadata()`](https://avisrilab.github.io/Matisse/reference/AddIsoformMetadata.md)
  : Add columns to the cell metadata

## Normalisation

Normalise count data before clustering. SCTransform applies variance
stabilisation and scales for sequencing depth. In event mode it
normalises the transcript assay; in junction mode it normalises the
gene-expression assay by default.

- [`SCTransform(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/SCTransform.MatisseObject.md)
  : SCTransform normalisation for MatisseObjects

## Calculate splicing ratios (PSI)

Compute PSI (Percent Spliced In) — the fraction of transcripts in each
cell that include a given exon. Values range from 0 (exon always
skipped) to 1 (exon always included). In event mode this is done at
construction; in junction mode call CalculatePSI() explicitly.

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

Plot splicing patterns across your cells. Overlay any feature — PSI
values, junction counts, or gene expression — on a UMAP, compare
splicing between cell types, or inspect junction usage for a gene of
interest. Pass the feature name via the `feature` argument.

- [`PlotUMAP()`](https://avisrilab.github.io/Matisse/reference/PlotUMAP.md)
  : UMAP plot coloured by any feature
- [`PlotViolin()`](https://avisrilab.github.io/Matisse/reference/PlotViolin.md)
  : Violin plot of feature values split by cell group
- [`PlotHeatmap()`](https://avisrilab.github.io/Matisse/reference/PlotHeatmap.md)
  : Heatmap of feature values across cells and events
- [`PlotCoverage()`](https://avisrilab.github.io/Matisse/reference/PlotCoverage.md)
  : Junction coverage bar plot for a gene
- [`CoveragePlot()`](https://avisrilab.github.io/Matisse/reference/CoveragePlot.md)
  : Sashimi-style coverage plot for a splice event
- [`PlotQCMetrics()`](https://avisrilab.github.io/Matisse/reference/PlotQCMetrics.md)
  : Violin/ridge plot of isoform QC metrics

## Utilities

Helper functions for building event tables and combining datasets.

- [`BuildSimpleEvents()`](https://avisrilab.github.io/Matisse/reference/BuildSimpleEvents.md)
  : Build a minimal junction event annotation table
- [`MergeMatisse()`](https://avisrilab.github.io/Matisse/reference/MergeMatisse.md)
  : Merge two MatisseObjects by cells
