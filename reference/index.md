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

Normalise count data before clustering. `SCTransform` applies variance
stabilisation and scales for sequencing depth — in event mode it targets
the `transcript` assay automatically; in junction mode it targets the
gene-expression assay. `NormalizeData`, `ScaleData`, and
`FindVariableFeatures` are log-normalisation alternatives.

- [`SCTransform(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/SCTransform.MatisseObject.md)
  : SCTransform normalisation for MatisseObjects
- [`NormalizeData(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/NormalizeData.MatisseObject.md)
  : Normalise gene-expression counts for a MatisseObject
- [`ScaleData(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/ScaleData.MatisseObject.md)
  : Scale gene-expression data for a MatisseObject
- [`FindVariableFeatures(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/FindVariableFeatures.MatisseObject.md)
  : Identify highly variable features for a MatisseObject

## Dimensionality reduction

Reduce the high-dimensional feature space to a compact representation
before clustering and visualisation. `RunPCA` is the standard route
after `SCTransform`; `RunSVD` (LSI) is used for ATAC-seq data in
multiome experiments.

- [`RunPCA(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/RunPCA.MatisseObject.md)
  : Run PCA on a MatisseObject
- [`RunUMAP(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/RunUMAP.MatisseObject.md)
  : Run UMAP on a MatisseObject
- [`RunTSNE(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/RunTSNE.MatisseObject.md)
  : Run t-SNE on a MatisseObject

## Clustering and differential expression

Group cells by transcriptional state and identify marker genes. All
functions operate on the embedded Seurat object and return the updated
`MatisseObject` — except `FindMarkers`, which returns a data frame.

- [`FindNeighbors(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/FindNeighbors.MatisseObject.md)
  : Compute a shared nearest-neighbour graph for a MatisseObject
- [`FindClusters(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/FindClusters.MatisseObject.md)
  : Cluster cells in a MatisseObject
- [`FindMarkers(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/FindMarkers.MatisseObject.md)
  : Find differentially expressed markers for a MatisseObject

## Metadata

Add or update per-cell metadata columns. New columns are accessible
immediately via
[`MatisseMeta()`](https://avisrilab.github.io/Matisse/reference/MatisseMeta.md)
and the `$` operator.

- [`AddMetaData(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/AddMetaData.MatisseObject.md)
  : Add metadata columns to a MatisseObject

## Signac / ATAC-seq methods

For multiome (paired ATAC + RNA) datasets. `RunTFIDF` normalises peak
counts; `FindTopFeatures` selects variable peaks; `RunSVD` performs
Latent Semantic Indexing (LSI) — the ATAC equivalent of PCA.

- [`RunTFIDF(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/RunTFIDF.MatisseObject.md)
  : Run TF-IDF normalisation for a MatisseObject
- [`FindTopFeatures(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/FindTopFeatures.MatisseObject.md)
  : Find highly variable ATAC-seq features for a MatisseObject
- [`RunSVD(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/RunSVD.MatisseObject.md)
  : Run SVD (LSI) for a MatisseObject

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
- [`PlotSashimi()`](https://avisrilab.github.io/Matisse/reference/PlotSashimi.md)
  : Sashimi-style coverage plot for a splice event
- [`PlotQCMetrics()`](https://avisrilab.github.io/Matisse/reference/PlotQCMetrics.md)
  : Violin/ridge plot of isoform QC metrics

## Utilities

Helper functions for building event tables and combining datasets.

- [`BuildSimpleEvents()`](https://avisrilab.github.io/Matisse/reference/BuildSimpleEvents.md)
  : Build a minimal junction event annotation table
- [`MergeMatisse()`](https://avisrilab.github.io/Matisse/reference/MergeMatisse.md)
  : Merge two MatisseObjects by cells
