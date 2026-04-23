# Identify highly variable features for a MatisseObject

Runs
[`FindVariableFeatures`](https://satijalab.org/seurat/reference/FindVariableFeatures.html)
on the embedded Seurat object and returns the updated `MatisseObject`.
Identifies genes whose expression varies most across cells — used to
select features for PCA.

## Usage

``` r
# S3 method for class 'MatisseObject'
FindVariableFeatures(object, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments forwarded to
  [`FindVariableFeatures`](https://satijalab.org/seurat/reference/FindVariableFeatures.html).

## Value

The updated `MatisseObject`.

## See also

[`RunPCA.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunPCA.MatisseObject.md),
[`NormalizeData.MatisseObject`](https://avisrilab.github.io/Matisse/reference/NormalizeData.MatisseObject.md)
