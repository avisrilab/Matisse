# Run SVD (LSI) for a MatisseObject

Runs [`RunSVD`](https://stuartlab.org/signac/reference/RunSVD.html) on
the embedded Seurat object and returns the updated `MatisseObject`.
Performs Latent Semantic Indexing (LSI) — the standard
dimensionality-reduction step for ATAC-seq data. Typically called after
[`RunTFIDF.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunTFIDF.MatisseObject.md).

## Usage

``` r
# S3 method for class 'MatisseObject'
RunSVD(object, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments forwarded to
  [`RunSVD`](https://stuartlab.org/signac/reference/RunSVD.html) (e.g.
  `n`, `assay`).

## Value

The updated `MatisseObject` with an `"lsi"` reduction.

## See also

[`RunTFIDF.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunTFIDF.MatisseObject.md),
[`RunUMAP.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunUMAP.MatisseObject.md)
