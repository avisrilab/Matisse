# Violin/ridge plot of isoform QC metrics

Violin/ridge plot of isoform QC metrics

Violin plots of isoform QC metrics

## Usage

``` r
PlotQCMetrics(object, features = NULL, ...)

# S4 method for class 'MatisseObject'
PlotQCMetrics(object, features = NULL, group_by = NULL, ncol = 2L)
```

## Arguments

- object:

  A `MatisseObject` with populated `isoform_metadata`.

- features:

  Character vector of QC metric names. Must be columns in
  `isoform_metadata`. Default: all numeric columns.

- ...:

  Additional arguments (see `PlotQCMetrics`).

- group_by:

  Character. Column in `Seurat::meta.data` to split cells by. Default:
  `NULL` (single group).

- ncol:

  Integer. Number of columns in the faceted output. Default: `2`.

## Value

A `ggplot` object.

A `ggplot` object.
