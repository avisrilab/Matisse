# Heatmap of feature values across cells and events

Heatmap of feature values across cells and events

Heatmap of PSI values (cells x events)

## Usage

``` r
PlotHeatmap(object, ...)

# S4 method for class 'MatisseObject'
PlotHeatmap(
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

  A `MatisseObject` with a `"psi"` assay.

- ...:

  Additional arguments (see `PlotHeatmap`).

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

A `ggplot` object.

A `ggplot` object.
