# Set the PSI matrix

Replaces the `"data"` layer of the `"psi"` `ChromatinAssay` inside the
embedded Seurat object.

## Usage

``` r
SetPSI(object, value)

# S4 method for class 'MatisseObject'
SetPSI(object, value)
```

## Arguments

- object:

  A `MatisseObject`.

- value:

  A sparse matrix (cells × events) of PSI values.

## Value

The updated `MatisseObject`.
