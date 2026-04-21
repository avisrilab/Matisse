# Summarize PSI distribution across cells for each event

Returns a summary table with per-event PSI statistics across all (or a
subset of) cells. Call this after
[`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
or after creating an event-mode object with
[`CreateMatisseObject`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md).

## Usage

``` r
SummarizePSI(object, cells = NULL)
```

## Arguments

- object:

  A `MatisseObject` with a `"psi"` assay.

- cells:

  Optional character vector of cell barcodes to subset.

## Value

A `data.frame` with one row per event and columns: `event_id`,
`mean_psi`, `median_psi`, `sd_psi`, `n_cells_covered`.
