# Violin plot of feature values split by cell group

Violin plot of feature values split by cell group

Violin plot of feature values split by cell group

## Usage

``` r
PlotViolin(
  object,
  feature,
  group_by = "seurat_clusters",
  colours = NULL,
  add_points = FALSE,
  title = NULL,
  ...
)

# S4 method for class 'MatisseObject'
PlotViolin(
  object,
  feature,
  group_by = "seurat_clusters",
  colours = NULL,
  add_points = FALSE,
  title = NULL,
  ...
)
```

## Arguments

- object:

  A `MatisseObject`.

- feature:

  Character. Feature to plot (PSI event ID, junction ID, or gene name).

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

  Character. Plot title. Defaults to the feature name.

## Value

A `ggplot` object.

A `ggplot` object.
