# Heatmap of PSI values (events x cells, DoHeatmap-style)

Draws a tile plot with splice events on the y-axis and cells on the
x-axis. Events are clustered by PSI profile. When `group_by` is
supplied, cells are ordered by group and a colored bar is drawn at the
top of the heatmap showing each cell's group (Seurat `DoHeatmap` style).
Both legends – continuous PSI and discrete groups – appear on the right.

## Usage

``` r
PlotHeatmap(
  object,
  events = NULL,
  cells = NULL,
  group_by = NULL,
  max_cells = 500L,
  max_events = 50L,
  na_colour = "grey90",
  title = NULL,
  ...
)

# S4 method for class 'MatisseObject'
PlotHeatmap(
  object,
  events = NULL,
  cells = NULL,
  group_by = NULL,
  max_cells = 500L,
  max_events = 50L,
  na_colour = "grey90",
  title = NULL,
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

  Character. Column in `Seurat::meta.data` used to order cells and draw
  the colored top bar. Default: `NULL`.

- max_cells:

  Integer. Downsample to this many cells before plotting. Default:
  `500`.

- max_events:

  Integer. Cap on events to plot. When the candidate set exceeds this,
  the top-variance events are selected automatically. Default: `50`,
  which keeps event-name labels readable; raise it for a denser heatmap
  (and expect smaller / overlapping y-axis labels).

- na_colour:

  Character. Colour for `NA` entries. Default: `"grey90"`.

- title:

  Character. Plot title. Default: `"PSI Heatmap"`.

- ...:

  Ignored; included for S4 generic compatibility.

## Value

A `ggplot` object.

A `ggplot` object.
