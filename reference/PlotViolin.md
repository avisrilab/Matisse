# Violin plot of feature values split by cell group

When `feature` is `NULL` (the default), automatically plots the QC
metrics `nCount_isoform`, `nFeature_isoform`, and `nPercent_isoform` as
a faceted panel.

Plots any feature – PSI event ID, junction ID, gene name, or metadata
column – as a violin/box split by a cell-level grouping variable. When
`feature = NULL` (the default), automatically plots the standard QC
metrics `nCount_isoform`, `nFeature_isoform`, and `nPercent_isoform` as
a faceted panel. Pass a character vector to `feature` to plot multiple
features in facets with free y-axis scales.

## Usage

``` r
PlotViolin(
  object,
  feature = NULL,
  group_by = "seurat_clusters",
  colours = NULL,
  add_points = FALSE,
  title = NULL,
  ncol = 2L,
  ...
)

# S4 method for class 'MatisseObject'
PlotViolin(
  object,
  feature = NULL,
  group_by = "seurat_clusters",
  colours = NULL,
  add_points = FALSE,
  title = NULL,
  ncol = 2L,
  ...
)
```

## Arguments

- object:

  A `MatisseObject`.

- feature:

  Character scalar or vector of features to plot (PSI event IDs,
  junction IDs, gene names, or metadata column names). `NULL` (the
  default) selects the standard QC metrics automatically.

- group_by:

  Character. Column in `Seurat::meta.data` to split cells by. Default:
  `"seurat_clusters"`.

- colours:

  Named character vector mapping group levels to colours. Default:
  `NULL` (uses ggplot2 defaults).

- add_points:

  Logical. Overlay individual cell values as jittered points. Default:
  `FALSE`.

- title:

  Character. Plot title. Defaults to the feature name when a single
  feature is requested; `NULL` for multiple features.

- ncol:

  Integer. Number of facet columns when `feature` is a vector of length
  \> 1. Default: `2`.

- ...:

  Ignored; included for S4 generic compatibility.

## Value

A `ggplot` object.

A `ggplot` object.
