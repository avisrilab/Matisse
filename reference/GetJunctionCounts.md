# Get raw junction count matrix (junction mode only)

Retrieves the per-junction read counts from the `"isoform"` `Assay5` in
junction mode.

## Usage

``` r
GetJunctionCounts(object, ...)

# S4 method for class 'MatisseObject'
GetJunctionCounts(object, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments (unused).

## Value

A sparse matrix (cells x junctions) of read counts, or `NULL` if the
object is in transcript mode.
