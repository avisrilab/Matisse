# Find highly variable ATAC-seq features for a MatisseObject

Runs
[`FindTopFeatures`](https://stuartlab.org/signac/reference/FindTopFeatures.html)
on the embedded Seurat object and returns the updated `MatisseObject`.
Selects the most accessible peaks for downstream LSI / UMAP.

## Usage

``` r
# S3 method for class 'MatisseObject'
FindTopFeatures(object, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments forwarded to
  [`FindTopFeatures`](https://stuartlab.org/signac/reference/FindTopFeatures.html)
  (e.g. `min.cutoff`).

## Value

The updated `MatisseObject` with top ATAC features flagged.

## See also

[`RunTFIDF.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunTFIDF.MatisseObject.md),
[`RunSVD.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunSVD.MatisseObject.md)
