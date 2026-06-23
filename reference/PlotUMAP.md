# UMAP plot – by group (DimPlot-style) or by feature (FeaturePlot-style)

Two modes, switched by `feature`:

- **Group mode** (default, `feature = NULL`): cells coloured discretely
  by a metadata column – like Seurat's `DimPlot`.

- **Feature mode** (`feature` supplied): cells coloured on a continuous
  scale by a PSI event ID, junction or transcript ID, gene name, or
  numeric metadata column – like Seurat's `FeaturePlot`. The palette
  adapts to the feature type (diverging RdBu \[0,1\] for PSI; sequential
  for counts and expression).

Two modes, switched by `feature`:

- **Group mode** (default, `feature = NULL`): cells coloured discretely
  by a metadata column – like Seurat's `DimPlot`. Pass `label = TRUE` to
  add group labels at cluster centroids.

- **Feature mode** (`feature` supplied): cells coloured on a continuous
  scale by a PSI event ID, junction or transcript ID, gene name, or
  numeric metadata column – like Seurat's `FeaturePlot`. The palette
  adapts to the feature type: diverging RdBu \[0,1\] for PSI; sequential
  for counts and expression.

## Usage

``` r
PlotUMAP(
  object,
  feature = NULL,
  reduction = "umap",
  dims = c(1L, 2L),
  group_by = "seurat_clusters",
  pt_size = 0.5,
  na_colour = "grey80",
  label = FALSE,
  title = NULL,
  ...
)

# S4 method for class 'MatisseObject'
PlotUMAP(
  object,
  feature = NULL,
  reduction = "umap",
  dims = c(1L, 2L),
  group_by = "seurat_clusters",
  pt_size = 0.5,
  na_colour = "grey80",
  label = FALSE,
  title = NULL,
  ...
)
```

## Arguments

- object:

  A `MatisseObject` with a UMAP reduction.

- feature:

  Character. Optional feature to plot – a PSI event ID (e.g.
  `"SE:chr1:100-200:300-400:+"`), a junction or transcript ID, a gene
  name, or a numeric cell-metadata column. `NULL` (the default) selects
  group mode.

- reduction:

  Character. Name of the dimensionality reduction to use. Default:
  `"umap"`.

- dims:

  Integer vector of length 2 selecting which dimensions to plot.
  Default: `c(1, 2)`.

- group_by:

  Character. Metadata column to colour by in group mode. Default:
  `"seurat_clusters"`. Ignored in feature mode.

- pt_size:

  Numeric. Point size. Default: `0.5`.

- na_colour:

  Character. Colour for cells with no data. Default: `"grey80"`.

- label:

  Logical. Group-mode only – add group labels at cluster centroids.
  Default: `FALSE`.

- title:

  Character. Plot title. Defaults to the feature name in feature mode
  and `NULL` in group mode.

- ...:

  Additional arguments (see `PlotUMAP`).

## Value

A `ggplot` object.

A `ggplot` object.
