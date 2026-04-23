# Normalise gene-expression counts for a MatisseObject

Runs
[`NormalizeData`](https://satijalab.org/seurat/reference/NormalizeData.html)
on the embedded Seurat object and returns the updated `MatisseObject`.
All splicing assays (junction, psi, transcript) are unaffected.

## Usage

``` r
# S3 method for class 'MatisseObject'
NormalizeData(object, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments forwarded to
  [`NormalizeData`](https://satijalab.org/seurat/reference/NormalizeData.html).

## Value

The updated `MatisseObject`.

## See also

[`SCTransform.MatisseObject`](https://avisrilab.github.io/Matisse/reference/SCTransform.MatisseObject.md),
[`ScaleData.MatisseObject`](https://avisrilab.github.io/Matisse/reference/ScaleData.MatisseObject.md)
