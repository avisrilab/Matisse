# Cluster cells in a MatisseObject

Runs
[`FindClusters`](https://satijalab.org/seurat/reference/FindClusters.html)
on the embedded Seurat object and returns the updated `MatisseObject`.
Cluster assignments are stored in `MatisseMeta(obj)$seurat_clusters` and
are immediately available for
[`PlotUMAP`](https://avisrilab.github.io/Matisse/reference/PlotUMAP.md)
and
[`PlotViolin`](https://avisrilab.github.io/Matisse/reference/PlotViolin.md).

## Usage

``` r
# S3 method for class 'MatisseObject'
FindClusters(object, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments forwarded to
  [`FindClusters`](https://satijalab.org/seurat/reference/FindClusters.html)
  (e.g. `resolution`).

## Value

The updated `MatisseObject` with `seurat_clusters` added to cell
metadata.

## See also

[`FindNeighbors.MatisseObject`](https://avisrilab.github.io/Matisse/reference/FindNeighbors.MatisseObject.md),
[`PlotUMAP`](https://avisrilab.github.io/Matisse/reference/PlotUMAP.md)
