# Merge two MatisseObjects by cells

Concatenates two `MatisseObject`s that share the same set of splice
events. All matrices are row-bound; the embedded Seurat objects are
merged via [`merge()`](https://rdrr.io/r/base/merge.html) (dispatches to
Seurat's merge method).

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
