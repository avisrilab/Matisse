# Run PCA on a MatisseObject

Runs [`RunPCA`](https://satijalab.org/seurat/reference/RunPCA.html) on
the embedded Seurat object and returns the updated `MatisseObject`. The
PCA result is stored inside the Seurat object and accessible via
`GetSeurat(obj)`.

## Usage

``` r
# S3 method for class 'MatisseObject'
RunPCA(object, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments forwarded to
  [`RunPCA`](https://satijalab.org/seurat/reference/RunPCA.html) (e.g.
  `assay`, `npcs`, `features`).

## Value

The updated `MatisseObject` with a `"pca"` reduction.

## Details

Typical usage after
[`SCTransform.MatisseObject`](https://avisrilab.github.io/Matisse/reference/SCTransform.MatisseObject.md):

    obj <- RunPCA(obj, assay = "SCT", npcs = 50)

## See also

[`RunUMAP.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunUMAP.MatisseObject.md),
[`SCTransform.MatisseObject`](https://avisrilab.github.io/Matisse/reference/SCTransform.MatisseObject.md)
