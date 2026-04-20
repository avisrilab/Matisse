# Merge two MatisseObjects by cells

Concatenates two `MatisseObject`s that share the same set of splice
events. The embedded Seurat objects are merged via
[`merge()`](https://rdrr.io/r/base/merge.html) (which dispatches to
Seurat's merge method, combining all assays including `"psi"` and
`"transcript"`).

## Usage

``` r
MergeMatisse(x, y, add_cell_ids = c("x", "y"), verbose = TRUE)
```

## Arguments

- x:

  A `MatisseObject`.

- y:

  A `MatisseObject`.

- add_cell_ids:

  Character vector of length 2. Prefixes appended to cell barcodes to
  avoid collisions. Default: `c("x", "y")`.

- verbose:

  Logical. Default: `TRUE`.

## Value

A merged `MatisseObject`.
