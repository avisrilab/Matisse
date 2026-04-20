# Summarize PSI distribution across cells for each event

Summarize PSI distribution across cells for each event

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
