# Package index

## Creating a Matisse object

Start here. This function combines your Seurat object with splicing data
into a single Matisse object (and computes PSI) in one call. Pass
`junction_counts` for short-read (junction) mode, or `transcript_counts`
for long-read (transcript) mode, plus splice events via `events`
(data.frame or path to a SUPPA2 `.ioe` file).

- [`show(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
  [`dim(`*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
  [`` `[`( ``*`<MatisseObject>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
  [`` `[[`( ``*`<MatisseObject>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
  [`` `$`( ``*`<MatisseObject>`*`)`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
  : The MatisseObject S4 class
- [`CreateMatisseObject()`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md)
  : Create a MatisseObject

## Retrieve or update your data

Matisse-specific accessors. For raw count and annotation matrices, use
the native Seurat API on the embedded Seurat object:
e.g. `SeuratObject::GetAssayData(GetSeurat(obj)[[“isoform”]], “counts”)`
for junction or transcript counts, or `GetSeurat(obj)[[“psi”]][[]]` for
the per-event annotation table.

- [`GetSeurat()`](https://avisrilab.github.io/Matisse/reference/GetSeurat.md)
  : Get the embedded Seurat object
- [`GetPSI()`](https://avisrilab.github.io/Matisse/reference/GetPSI.md)
  : Get the PSI matrix
- [`SetPSI()`](https://avisrilab.github.io/Matisse/reference/SetPSI.md)
  : Set the PSI matrix
- [`MatisseMeta()`](https://avisrilab.github.io/Matisse/reference/MatisseMeta.md)
  [`` `MatisseMeta<-`() ``](https://avisrilab.github.io/Matisse/reference/MatisseMeta.md)
  : Get or set cell-level metadata

## Normalisation

Normalise count data before clustering. `SCTransform` applies variance
stabilisation and scales for sequencing depth – in transcript mode it
targets the `isoform` assay automatically; in junction mode it targets
the gene-expression assay. `NormalizeData`, `ScaleData`, and
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

Compute PSI (Percent Spliced In) – the fraction of transcripts in each
cell that include a given exon. Values range from 0 (exon always
skipped) to 1 (exon always included).
[`CreateMatisseObject()`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md)
calls
[`CalculatePSI()`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
automatically; you only call it directly to recompute with different
parameters (e.g. `min_coverage`) or after constructing with
`defer_psi = TRUE`.

- [`CalculatePSI()`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
  : Calculate PSI matrix from junction or transcript counts
- [`SummarizePSI()`](https://avisrilab.github.io/Matisse/reference/SummarizePSI.md)
  : Summarize PSI distribution across cells for each event

## Quality control and filtering

Identify and remove cells or splicing events that don’t have enough data
for reliable analysis. QC columns written automatically:
`nCount_isoform` and `nFeature_isoform` at construction;
`nPercent_isoform` after CalculatePSI().

- [`FilterCells()`](https://avisrilab.github.io/Matisse/reference/FilterCells.md)
  : Filter cells by QC thresholds
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
  : Heatmap of PSI values (events x cells, DoHeatmap style)
- [`PlotSashimi()`](https://avisrilab.github.io/Matisse/reference/PlotSashimi.md)
  : Sashimi-style coverage plot for a splice event

## Utilities

Helper functions for combining datasets.

- [`MergeMatisse()`](https://avisrilab.github.io/Matisse/reference/MergeMatisse.md)
  : Merge two MatisseObjects by cells
