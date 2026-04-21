# Junction coverage bar plot for a gene

Aggregates junction read counts across selected cells and plots a bar
chart ordered by genomic position. Only available for objects in
junction mode.

## Usage

``` r
PlotCoverage(object, gene, ...)

# S4 method for class 'MatisseObject'
PlotCoverage(object, gene, cells = NULL, log_scale = FALSE)
```

## Arguments

- object:

  A `MatisseObject` in junction mode with `junction_data` populated.

- gene:

  Character. Gene name to plot.

- ...:

  Additional arguments (see `PlotCoverage`).

- cells:

  Character vector of cell barcodes to aggregate over. Default: all
  cells.

- log_scale:

  Logical. Use log10 y-axis. Default: `FALSE`.

## Value

A `ggplot` object.

A `ggplot` object.
