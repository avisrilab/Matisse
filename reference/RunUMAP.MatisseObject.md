# Run UMAP on a MatisseObject

Runs [`RunUMAP`](https://satijalab.org/seurat/reference/RunUMAP.html) on
the embedded Seurat object and returns the updated `MatisseObject`. Call
after
[`RunPCA.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunPCA.MatisseObject.md)
(or
[`RunSVD.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunSVD.MatisseObject.md)
for ATAC data). The resulting embedding is accessible via
`GetSeurat(obj)` and used by
[`PlotUMAP`](https://avisrilab.github.io/Matisse/reference/PlotUMAP.md).

## Usage

``` r
# S3 method for class 'MatisseObject'
RunUMAP(object, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments forwarded to
  [`RunUMAP`](https://satijalab.org/seurat/reference/RunUMAP.html) (e.g.
  `dims`, `reduction`).

## Value

The updated `MatisseObject` with a `"umap"` reduction.

## See also

[`RunPCA.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunPCA.MatisseObject.md),
[`PlotUMAP`](https://avisrilab.github.io/Matisse/reference/PlotUMAP.md)
