# Violin/ridge plot of isoform QC metrics

Violin/ridge plot of isoform QC metrics

Violin plots of isoform QC metrics

## Usage

``` r
PlotQCMetrics(object, features = NULL, group_by = NULL, ncol = 2L, ...)

# S4 method for class 'MatisseObject'
PlotQCMetrics(object, features = NULL, group_by = NULL, ncol = 2L, ...)
```

## Arguments

- object:

  A `MatisseObject` with QC metrics computed by
  [`ComputeIsoformQC`](https://avisrilab.github.io/Matisse/reference/ComputeIsoformQC.md).

- features:

  Character vector of QC metric names. Must be columns in
  `MatisseMeta(object)`. Default: all numeric columns.

- group_by:

  Character. Column in `Seurat::meta.data` to split cells by. Default:
  `NULL` (single group).

- ncol:

  Integer. Number of columns in the faceted output. Default: `2`.

- ...:

  Additional arguments passed to methods.

## Value

A `ggplot` object.

A `ggplot` object.
