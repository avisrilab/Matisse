# Heatmap of feature values across cells and events

Heatmap of feature values across cells and events

Heatmap of PSI values (cells x events)

## Usage

``` r
PlotHeatmap(
  object,
  events = NULL,
  cells = NULL,
  group_by = NULL,
  max_cells = 500L,
  max_events = 200L,
  na_colour = "grey90",
  ...
)

# S4 method for class 'MatisseObject'
PlotHeatmap(
  object,
  events = NULL,
  cells = NULL,
  group_by = NULL,
  max_cells = 500L,
  max_events = 200L,
  na_colour = "grey90",
  ...
)
```

## Arguments

- object:

  A `MatisseObject` with a `"psi"` assay.

- events:

  Character vector of event IDs to include. Default: `NULL`
  (top-variance events up to `max_events`).

- cells:

  Character vector of cell barcodes to include. Default: `NULL` (random
  sample up to `max_cells`).

- group_by:

  Character. Column in `Seurat::meta.data` used to annotate and order
  cells. Default: `NULL`.

- max_cells:

  Integer. Downsample to this many cells before plotting. Default:
  `500`.

- max_events:

  Integer. Cap on events to plot. When the candidate set exceeds this,
  the top-variance events are selected automatically. Default: `200`.

- na_colour:

  Character. Colour for `NA` entries. Default: `"grey90"`.

## Value

A `ggplot` object.

A `ggplot` object.
