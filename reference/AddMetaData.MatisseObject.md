# Add metadata columns to a MatisseObject

Runs
[`AddMetaData`](https://satijalab.github.io/seurat-object/reference/AddMetaData.html)
on the embedded Seurat object and returns the updated `MatisseObject`.
New columns are immediately accessible via
[`MatisseMeta`](https://avisrilab.github.io/Matisse/reference/MatisseMeta.md)
and the `$` operator.

## Usage

``` r
# S3 method for class 'MatisseObject'
AddMetaData(object, metadata, col.name = NULL, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- metadata:

  A named vector, `data.frame`, or list of columns to add.

- col.name:

  Character. Column name to use when `metadata` is a vector. Default:
  `NULL`.

- ...:

  Additional arguments forwarded to
  [`AddMetaData`](https://satijalab.github.io/seurat-object/reference/AddMetaData.html).

## Value

The updated `MatisseObject` with new metadata columns.

## See also

[`MatisseMeta`](https://avisrilab.github.io/Matisse/reference/MatisseMeta.md),
[`AddIsoformMetadata`](https://avisrilab.github.io/Matisse/reference/AddIsoformMetadata.md)
