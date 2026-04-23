# Run TF-IDF normalisation for a MatisseObject

Runs [`RunTFIDF`](https://stuartlab.org/signac/reference/RunTFIDF.html)
on the embedded Seurat object and returns the updated `MatisseObject`.
Used for ATAC-seq peak counts in multiome datasets before
[`RunSVD.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunSVD.MatisseObject.md).

## Usage

``` r
# S3 method for class 'MatisseObject'
RunTFIDF(object, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments forwarded to
  [`RunTFIDF`](https://stuartlab.org/signac/reference/RunTFIDF.html)
  (e.g. `assay`, `method`).

## Value

The updated `MatisseObject` with TF-IDF normalised counts.

## See also

[`RunSVD.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunSVD.MatisseObject.md),
[`FindTopFeatures.MatisseObject`](https://avisrilab.github.io/Matisse/reference/FindTopFeatures.MatisseObject.md)
