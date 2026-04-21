# Get or set cell-level metadata

Returns the full `meta.data` of the embedded Seurat object, which
includes all per-cell QC metrics and annotations added by Matisse (e.g.
`n_junctions_detected`, `mean_psi`) alongside standard Seurat columns.
Use
[`AddIsoformMetadata()`](https://avisrilab.github.io/Matisse/reference/AddIsoformMetadata.md)
to add new columns.

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
