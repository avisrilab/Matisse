# Find highly variable ATAC-seq features for a MatisseObject

Runs
[`FindTopFeatures`](https://stuartlab.org/signac/reference/FindTopFeatures.html)
on the embedded Seurat object and returns the updated `MatisseObject`.
Selects the most accessible peaks for downstream LSI / UMAP.

## Usage

``` r
# S3 method for class 'MatisseObject'
FindTopFeatures(object, assay = NULL, min.cutoff = "q5", verbose = TRUE, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- assay:

  Name of assay to use

- min.cutoff:

  Cutoff for feature to be included in the VariableFeatures for the
  object. This can be a percentile specified as 'q' followed by the
  minimum percentile, for example 'q5' to set the top 95\\ as the
  VariableFeatures for the object. Alternatively, this can be an integer
  specifying the minimum number of counts for the feature to be included
  in the set of VariableFeatures. For example, setting to 10 will
  include features with \>10 total counts in the set of
  VariableFeatures. If NULL, include all features in VariableFeatures.
  If NA, VariableFeatures will not be altered, and only the feature
  metadata will be updated with the total counts and percentile rank for
  each feature.

- verbose:

  Display messages

- ...:

  Additional arguments forwarded to
  [`FindTopFeatures`](https://stuartlab.org/signac/reference/FindTopFeatures.html)
  (e.g. `min.cutoff`).

## Value

The updated `MatisseObject` with top ATAC features flagged.

## See also

[`RunTFIDF.MatisseObject`](https://avisrilab.org/Matisse/reference/RunTFIDF.MatisseObject.md),
[`RunSVD.MatisseObject`](https://avisrilab.org/Matisse/reference/RunSVD.MatisseObject.md)
