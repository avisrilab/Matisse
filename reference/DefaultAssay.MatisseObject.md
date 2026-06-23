# Get or set the default assay of a MatisseObject

Delegates to
[`DefaultAssay`](https://satijalab.github.io/seurat-object/reference/DefaultAssay.html)
on the embedded Seurat object. Lets users switch between `"RNA"`,
`"isoform"`, `"psi"`, `"SCT"`, etc. without unwrapping the object first.

## Usage

``` r
# S3 method for class 'MatisseObject'
DefaultAssay(object, ...)

# S3 method for class 'MatisseObject'
DefaultAssay(object, ...) <- value
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments forwarded to
  [`DefaultAssay`](https://satijalab.github.io/seurat-object/reference/DefaultAssay.html).

- value:

  Character. The new default assay name (setter only).

## Value

The default assay name (getter) or the updated `MatisseObject` (setter).
