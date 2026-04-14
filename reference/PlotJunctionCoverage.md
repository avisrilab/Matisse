# Per-cell junction coverage plot for a gene

Displays total junction read counts across all cells for each junction
in the specified gene, grouped by event type.

## Usage

``` r
PlotJunctionCoverage(object, gene, ...)

# S4 method for class 'MatisseObject'
PlotJunctionCoverage(object, gene, cells = NULL, log_scale = FALSE)
```

## Arguments

- object:

  A `MatisseObject` with non-`NULL` `junction_counts` and
  `junction_data` slots.

- gene:

  Character. Gene name to plot.

- ...:

  Additional arguments (see `PlotJunctionCoverage`).

- cells:

  Character vector of cell barcodes to aggregate over. Default: all
  cells.

- log_scale:

  Logical. Use log10 y-axis. Default: `FALSE`.

## Value

A `ggplot` object.

A `ggplot` object.
