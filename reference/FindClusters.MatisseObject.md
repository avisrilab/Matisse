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
FindClusters(
  object,
  graph.name = NULL,
  cluster.name = NULL,
  modularity.fxn = 1,
  initial.membership = NULL,
  node.sizes = NULL,
  resolution = 0.8,
  method = NULL,
  algorithm = 1,
  leiden_method = c("leidenbase", "igraph"),
  leiden_objective_function = c("modularity", "CPM"),
  n.start = 10,
  n.iter = 10,
  random.seed = 0,
  group.singletons = TRUE,
  temp.file.location = NULL,
  edge.file.name = NULL,
  verbose = TRUE,
  ...
)
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
