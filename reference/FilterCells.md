# Filter cells by QC thresholds

Removes cells that do not pass the specified thresholds on QC columns in
`MatisseMeta(object)` (the Seurat `meta.data`). The columns
`nCount_isoform` and `nFeature_isoform` are written at construction by
[`CreateMatisseObject`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md);
`nPercent_isoform` is written by
[`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md).

## Usage

``` r
FilterCells(
  object,
  min_features_isoform = NULL,
  max_features_isoform = NULL,
  min_counts_isoform = NULL,
  max_counts_isoform = NULL,
  min_pct_isoform = NULL,
  custom_filters = NULL,
  verbose = TRUE,
  ...
)

# S4 method for class 'MatisseObject'
FilterCells(
  object,
  min_features_isoform = NULL,
  max_features_isoform = NULL,
  min_counts_isoform = NULL,
  max_counts_isoform = NULL,
  min_pct_isoform = NULL,
  custom_filters = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A `MatisseObject`.

- min_features_isoform:

  Integer. Minimum `nFeature_isoform` (number of junctions or
  transcripts with at least one read). Default: `NULL`.

- max_features_isoform:

  Integer. Maximum `nFeature_isoform`. Default: `NULL`.

- min_counts_isoform:

  Integer. Minimum `nCount_isoform` (total reads in the isoform assay).
  Default: `NULL`.

- max_counts_isoform:

  Integer. Maximum `nCount_isoform`. Default: `NULL`.

- min_pct_isoform:

  Numeric (0-100). Minimum `nPercent_isoform` (percentage of splice
  events with a non-NA PSI value). Requires
  [`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
  to have been run. Default: `NULL`.

- custom_filters:

  Named list of two-element numeric vectors `c(min, max)` applied to
  arbitrary metadata columns. Use `NA` for a one-sided bound. Default:
  `NULL`.

- verbose:

  Logical. Default: `TRUE`.

- ...:

  Ignored; included for S4 generic compatibility.

## Value

The filtered `MatisseObject`.

The filtered `MatisseObject`.

## See also

[`FilterEvents`](https://avisrilab.github.io/Matisse/reference/FilterEvents.md),
[`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
