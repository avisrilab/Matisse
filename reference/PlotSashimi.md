# Sashimi-style coverage plot for a splice event

Draws junction arcs scaled by read count over a schematic gene
structure. Arcs are coloured by role: inclusion (blue) vs exclusion
(red). Works in both junction mode and transcript mode. Optionally
faceted by a cell metadata column.

Draws junction arcs scaled by aggregate read count over a schematic gene
structure. Arcs are coloured by role: inclusion (blue) vs exclusion
(red).

## Usage

``` r
PlotSashimi(
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
PlotSashimi(
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

  Character. A single event ID. Run `rownames(GetSeurat(obj)[["psi"]])`
  to list available IDs; typical SE format is
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

  Additional arguments (see `PlotSashimi`).

## Value

A `ggplot` object.

A `ggplot` object.

## Details

In **junction mode** each arc corresponds to an individual junction with
its own read count. In **transcript mode** the SE event_id is parsed to
derive junction coordinates; inclusion and exclusion counts come from
the `"counts"` and `"exclusion"` layers of the PSI assay.

Supported event types in transcript mode: **SE** (skipped exon) and
**RI** (retained intron). Junction mode supports all event types since
coordinates are auto-parsed from junction IDs.

**Note on transcript-mode SE arcs:** transcript-level counting
aggregates reads to events, not to individual junctions. The total
inclusion-event count is therefore split evenly across the two SE
inclusion arcs in the plot. The two arcs may have differed in reality;
use junction-mode input if you need per-junction read counts.
