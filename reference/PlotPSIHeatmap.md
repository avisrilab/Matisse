# Heatmap of PSI values across cells and events

Renders a heatmap of PSI values. Cells are optionally ordered by a
metadata column; events are ordered by hierarchical clustering.

## Usage

``` r
PlotPSIHeatmap(object, ...)

# S4 method for class 'MatisseObject'
PlotPSIHeatmap(
  object,
  events = NULL,
  cells = NULL,
  group_by = NULL,
  max_cells = 500L,
  na_colour = "grey90"
)
```

## Arguments

- object:

  A `MatisseObject` with a non-`NULL` `psi` slot.

- ...:

  Additional arguments (see `PlotPSIHeatmap`).

- events:

  Character vector of event IDs to include. Default: all events.

- cells:

  Character vector of cell barcodes to include. Default: all cells.

- group_by:

  Character. Column in `Seurat::meta.data` used to annotate and order
  cells. Default: `NULL`.

- max_cells:

  Integer. Downsample to this many cells before plotting. Default:
  `500`.

- na_colour:

  Character. Colour for `NA` entries. Default: `"grey90"`.

## Value

A `ggplot` or `pheatmap` object.

A `ggplot` object.
