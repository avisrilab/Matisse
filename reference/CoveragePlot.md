# Sashimi-style coverage plot for a splice event

Draws junction arcs scaled by read count over a schematic gene
structure. Arcs are coloured by role: inclusion (blue) vs exclusion
(red). Works in both junction mode (per-junction counts) and event mode
(aggregated inclusion/exclusion counts). Optionally faceted by a cell
metadata column.

Draws junction arcs scaled by aggregate read count over a schematic gene
structure. Arcs are coloured by role: inclusion (blue) vs exclusion
(red).

## Usage

``` r
CoveragePlot(
  object,
  event_id,
  cells = NULL,
  group_by = NULL,
  arc_scale = c("sqrt", "linear", "log"),
  colours = c(inclusion = "#4393c3", exclusion = "#d6604d"),
  title = NULL,
  ...
)

# S4 method for class 'MatisseObject'
CoveragePlot(
  object,
  event_id,
  cells = NULL,
  group_by = NULL,
  arc_scale = c("sqrt", "linear", "log"),
  colours = c(inclusion = "#4393c3", exclusion = "#d6604d"),
  title = NULL,
  ...
)
```

## Arguments

- object:

  A `MatisseObject` with a PSI assay computed.

- event_id:

  Character. Event ID as stored in `event_data`, e.g.
  `"SE:chr1:1201-2999:3201-4999:+"`.

- cells:

  Character vector of cell barcodes to aggregate over. Default: all
  cells.

- group_by:

  Character. Column in Seurat meta.data to facet by. Default: `NULL`
  (all cells pooled).

- arc_scale:

  Character. How to scale arc height to read count: `"sqrt"` (default),
  `"linear"`, or `"log"`.

- colours:

  Named character vector with elements `"inclusion"` and `"exclusion"`
  giving arc colours.

- title:

  Character. Plot title. Defaults to `event_id`.

- ...:

  Additional arguments (see `CoveragePlot`).

## Value

A `ggplot` object.

A `ggplot` object.

## Details

In **junction mode** each arc corresponds to an individual junction with
its own read count. In **event mode** the SE event_id is parsed to
derive junction coordinates; inclusion and exclusion counts come from
the `"counts"` and `"exclusion"` layers of the PSI assay.

Currently only SE (skipped exon) events are supported for coordinate
derivation in event mode. Other event types require junction mode.
