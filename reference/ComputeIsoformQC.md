# Compute per-cell isoform QC metrics

Calculates a panel of QC metrics from the junction count and PSI layers
and stores them directly in the embedded Seurat object's `meta.data`.
Existing columns with the same names are overwritten.

## Usage

``` r
ComputeIsoformQC(object, ...)

# S4 method for class 'MatisseObject'
ComputeIsoformQC(object, min_coverage = 5L, verbose = TRUE)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments (unused).

- min_coverage:

  Integer. Minimum total reads to consider an event "covered" when
  computing `n_events_covered`. Default: `5`.

- verbose:

  Logical. Default: `TRUE`.

## Value

The updated `MatisseObject` with QC columns added to cell metadata.

The `MatisseObject` with QC columns added to `MatisseMeta(object)` (i.e.
the Seurat `meta.data`).

## Details

Computed metrics:

- n_junctions_detected:

  Number of junctions with at least one read (junction mode only).

- total_junction_reads:

  Total junction read count across all junctions (junction mode only).

- n_events_covered:

  Number of splice events with PSI \\\geq\\ `min_coverage` (requires PSI
  to be calculated first).

- pct_events_covered:

  Percentage of events with sufficient coverage.

- mean_psi:

  Mean PSI across covered events.

## See also

[`FilterCells`](https://avisrilab.github.io/Matisse/reference/FilterCells.md),
[`PlotQCMetrics`](https://avisrilab.github.io/Matisse/reference/PlotQCMetrics.md)
