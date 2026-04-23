# Scale gene-expression data for a MatisseObject

Runs
[`ScaleData`](https://satijalab.org/seurat/reference/ScaleData.html) on
the embedded Seurat object and returns the updated `MatisseObject`.
Typically called after
[`NormalizeData.MatisseObject`](https://avisrilab.github.io/Matisse/reference/NormalizeData.MatisseObject.md)
or
[`FindVariableFeatures.MatisseObject`](https://avisrilab.github.io/Matisse/reference/FindVariableFeatures.MatisseObject.md).

## Usage

``` r
# S3 method for class 'MatisseObject'
ScaleData(object, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments forwarded to
  [`ScaleData`](https://satijalab.org/seurat/reference/ScaleData.html).

## Value

The updated `MatisseObject`.

## See also

[`NormalizeData.MatisseObject`](https://avisrilab.github.io/Matisse/reference/NormalizeData.MatisseObject.md),
[`SCTransform.MatisseObject`](https://avisrilab.github.io/Matisse/reference/SCTransform.MatisseObject.md)
