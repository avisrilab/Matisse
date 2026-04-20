# Violin plot of PSI values split by cell group

Violin plot of PSI values split by cell group

Violin plot of PSI values split by cell group

## Usage

``` r
PlotPSIViolin(object, event_id, ...)

# S4 method for class 'MatisseObject'
PlotPSIViolin(
  object,
  event_id,
  group_by = "seurat_clusters",
  colours = NULL,
  add_points = FALSE,
  title = NULL
)
```

## Arguments

- object:

  A `MatisseObject` with a `"psi"` assay.

- event_id:

  Character. Column name in the PSI matrix.

- ...:

  Additional arguments (see `PlotPSIViolin`).

- group_by:

  Character. Column in `Seurat::meta.data` to split cells by. Default:
  `"seurat_clusters"`.

- colours:

  Named character vector mapping group levels to colours. Default:
  `NULL` (uses ggplot2 defaults).

- add_points:

  Logical. Overlay individual cell PSI values as jittered points.
  Default: `FALSE`.

- title:

  Character. Plot title. Defaults to the event ID.

## Value

A `ggplot` object.

A `ggplot` object.
