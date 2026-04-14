# Filter cells by isoform QC thresholds

Removes cells that do not pass the specified thresholds on columns in
the `isoform_metadata` slot. Pass named minimum and/or maximum bounds.

## Usage

``` r
FilterCells(object, ...)

# S4 method for class 'MatisseObject'
FilterCells(
  object,
  min_junctions = NULL,
  max_junctions = NULL,
  min_junction_reads = NULL,
  max_junction_reads = NULL,
  min_pct_covered = NULL,
  custom_filters = NULL,
  verbose = TRUE
)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Named numeric thresholds; see `FilterCells`.

- min_junctions:

  Integer. Minimum `n_junctions_detected`. Default: `NULL` (no filter).

- max_junctions:

  Integer. Maximum `n_junctions_detected`. Default: `NULL`.

- min_junction_reads:

  Integer. Minimum `total_junction_reads`. Default: `NULL`.

- max_junction_reads:

  Integer. Maximum `total_junction_reads`. Default: `NULL`.

- min_pct_covered:

  Numeric (0–100). Minimum `pct_events_covered`. Default: `NULL`.

- custom_filters:

  Named list of two-element numeric vectors `c(min, max)` applied to
  arbitrary `isoform_metadata` columns. Use `NA` for a one-sided bound,
  e.g. `list(mean_psi = c(0.1, NA))`. Default: `NULL`.

- verbose:

  Logical. Default: `TRUE`.

## Value

The filtered `MatisseObject`.

The filtered `MatisseObject`.

## See also

[`ComputeIsoformQC`](https://avisrilab.github.io/Matisse/reference/ComputeIsoformQC.md)
