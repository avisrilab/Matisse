# Compute a shared nearest-neighbour graph for a MatisseObject

Runs
[`FindNeighbors`](https://satijalab.org/seurat/reference/FindNeighbors.html)
on the embedded Seurat object and returns the updated `MatisseObject`.
Typically called after
[`RunPCA.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunPCA.MatisseObject.md)
and before
[`FindClusters.MatisseObject`](https://avisrilab.github.io/Matisse/reference/FindClusters.MatisseObject.md).

## Usage

``` r
# S3 method for class 'MatisseObject'
FindNeighbors(object, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments forwarded to
  [`FindNeighbors`](https://satijalab.org/seurat/reference/FindNeighbors.html)
  (e.g. `dims`, `reduction`).

## Value

The updated `MatisseObject` with neighbour graph stored inside the
embedded Seurat object.

## See also

[`FindClusters.MatisseObject`](https://avisrilab.github.io/Matisse/reference/FindClusters.MatisseObject.md),
[`RunPCA.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunPCA.MatisseObject.md)
