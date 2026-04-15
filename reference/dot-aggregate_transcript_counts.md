# Aggregate transcript counts into per-event inclusion/exclusion matrices

Aggregate transcript counts into per-event inclusion/exclusion matrices

## Usage

``` r
.aggregate_transcript_counts(tx_counts, events, min_coverage, cells)
```

## Arguments

- tx_counts:

  Matrix (transcripts x cells).

- events:

  data.frame from `.parse_ioe_files`.

- min_coverage:

  Integer coverage threshold.

- cells:

  Character vector of cell barcodes (column order).

## Value

List with elements `psi`, `inclusion`, `exclusion` (each a sparse cells
x events matrix).
