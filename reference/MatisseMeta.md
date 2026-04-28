# Get or set cell-level metadata

Returns the full `meta.data` of the embedded Seurat object, which
includes all per-cell QC metrics and annotations added by Matisse (e.g.
`nCount_isoform`, `nFeature_isoform`, `nPercent_isoform`) alongside
standard Seurat columns. To add new columns, use Seurat's
[`AddMetaData`](https://satijalab.github.io/seurat-object/reference/AddMetaData.html)
(an S3 method dispatches on `MatisseObject`).

## Usage

``` r
MatisseMeta(object, ...)

MatisseMeta(object) <- value

# S4 method for class 'MatisseObject'
MatisseMeta(object, ...)

# S4 method for class 'MatisseObject'
MatisseMeta(object) <- value
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments (unused).

- value:

  A `data.frame` whose columns are added to cell metadata (for the
  setter). Rownames must match cell barcodes.

## Value

For the getter: a `data.frame`. For the setter: the updated
`MatisseObject`.
