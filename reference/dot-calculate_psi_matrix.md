# Compute PSI matrix from a junction count matrix

Compute PSI matrix from a junction count matrix

## Usage

``` r
.calculate_psi_matrix(jxn_counts, events, min_coverage, na_fill, verbose)
```

## Arguments

- jxn_counts:

  Sparse matrix (cells x junctions).

- events:

  data.frame with event definitions.

- min_coverage:

  Integer threshold.

- na_fill:

  Fill value for low-coverage entries.

- verbose:

  Logical.

## Value

A list with elements `psi`, `inclusion`, `exclusion`.
