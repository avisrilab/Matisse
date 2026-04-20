# Get the PSI matrix

Retrieves the PSI (Percent Spliced In) matrix from the `"psi"`
`ChromatinAssay` stored inside the embedded Seurat object.

## Usage

``` r
GetPSI(object, ...)

# S4 method for class 'MatisseObject'
GetPSI(object, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments (unused).

## Value

A sparse matrix (cells × events) of PSI values in \\\[0,1\]\\. `NULL` if
no `"psi"` assay exists yet.
