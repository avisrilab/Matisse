# Find differentially expressed markers for a MatisseObject

Runs
[`FindMarkers`](https://satijalab.org/seurat/reference/FindMarkers.html)
on the embedded Seurat object. Unlike most dispatch methods, this
returns a `data.frame` of marker statistics rather than an updated
`MatisseObject`.

## Usage

``` r
# S3 method for class 'MatisseObject'
FindMarkers(object, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments forwarded to
  [`FindMarkers`](https://satijalab.org/seurat/reference/FindMarkers.html)
  (e.g. `ident.1`, `ident.2`, `group.by`, `features`).

## Value

A `data.frame` of marker genes with columns `p_val`, `avg_log2FC`,
`pct.1`, `pct.2`, `p_val_adj`.

## See also

[`FindClusters.MatisseObject`](https://avisrilab.github.io/Matisse/reference/FindClusters.MatisseObject.md)
