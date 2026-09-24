# Run TF-IDF normalisation for a MatisseObject

Runs [`RunTFIDF`](https://stuartlab.org/signac/reference/RunTFIDF.html)
on the embedded Seurat object and returns the updated `MatisseObject`.
Used for ATAC-seq peak counts in multiome datasets before
[`RunSVD.MatisseObject`](https://avisrilab.org/Matisse/reference/RunSVD.MatisseObject.md).

## Usage

``` r
# S3 method for class 'MatisseObject'
RunTFIDF(
  object,
  assay = NULL,
  method = 1,
  scale.factor = 10000,
  idf = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A `MatisseObject`.

- assay:

  Name of assay to use

- method:

  Which TF-IDF implementation to use. Choice of:

  - 1: The TF-IDF implementation used by Stuart & Butler et al. 2019
    ([doi:10.1101/460147](https://doi.org/10.1101/460147)). This
    computes \\\log(TF \times IDF)\\.

  - 2: The TF-IDF implementation used by Cusanovich & Hill et al. 2018
    ([doi:10.1016/j.cell.2018.06.052](https://doi.org/10.1016/j.cell.2018.06.052)).
    This computes \\TF \times (\log(IDF))\\.

  - 3: The log-TF method used by Andrew Hill. This computes \\\log(TF)
    \times \log(IDF)\\.

  - 4: The 10x Genomics method (no TF normalization). This computes
    \\IDF\\.

- scale.factor:

  Which scale factor to use. Default is 10000.

- idf:

  A precomputed IDF vector to use. If NULL, compute based on the input
  data matrix.

- verbose:

  Print progress

- ...:

  Additional arguments forwarded to
  [`RunTFIDF`](https://stuartlab.org/signac/reference/RunTFIDF.html)
  (e.g. `assay`, `method`).

## Value

The updated `MatisseObject` with TF-IDF normalised counts.

## See also

[`RunSVD.MatisseObject`](https://avisrilab.org/Matisse/reference/RunSVD.MatisseObject.md),
[`FindTopFeatures.MatisseObject`](https://avisrilab.org/Matisse/reference/FindTopFeatures.MatisseObject.md)
