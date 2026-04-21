# Get raw junction count matrix

Retrieves the per-junction read counts from the `"junction"` `Assay5`
stored inside the embedded Seurat object.

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

A sparse matrix (cells × junctions) of read counts, or `NULL` if the
object is in event mode or no junction assay exists.
