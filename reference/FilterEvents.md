# Filter splice events by coverage or variance

Removes events that do not pass minimum coverage or variance thresholds.
Operates on the PSI matrix.

## Usage

``` r
FilterEvents(object, ...)

# S4 method for class 'MatisseObject'
FilterEvents(
  object,
  min_cells_covered = 10L,
  min_psi_variance = NULL,
  verbose = TRUE
)
```

## Arguments

- object:

  A `MatisseObject` with a non-`NULL` `psi` slot.

- ...:

  Named thresholds; see `FilterEvents`.

- min_cells_covered:

  Integer. Minimum number of cells in which the event must have a
  non-`NA` PSI value. Default: `10`.

- min_psi_variance:

  Numeric. Minimum variance of PSI across covered cells. Filters out
  constitutively spliced (near PSI = 0 or 1) events. Default: `NULL` (no
  variance filter).

- verbose:

  Logical. Default: `TRUE`.

## Value

The filtered `MatisseObject`.

The filtered `MatisseObject`.
