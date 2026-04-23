# SCTransform normalisation for MatisseObjects

Runs
[`SCTransform`](https://satijalab.org/seurat/reference/SCTransform.html)
with mode-aware defaults. In **event mode**, normalises the
`"transcript"` assay. In **junction mode**, normalises the active
default assay (usually `"RNA"`). Override with the `assay` argument.

## Usage

``` r
# S3 method for class 'MatisseObject'
SCTransform(object, assay = NULL, vars_to_regress = NULL, verbose = TRUE, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- assay:

  Character. Assay to normalise. Default: `"transcript"` in event mode;
  the active default assay in junction mode.

- vars_to_regress:

  Character vector. Variables to regress out. Default: `NULL`.

- verbose:

  Logical. Default: `TRUE`.

- ...:

  Additional arguments forwarded to
  [`SCTransform`](https://satijalab.org/seurat/reference/SCTransform.html).

## Value

The updated `MatisseObject`.

## Details

PCA is **not** run automatically. Call
[`RunPCA()`](https://satijalab.org/seurat/reference/RunPCA.html) on the
resulting object after normalisation:

    obj <- SCTransform(obj)
    obj <- RunPCA(obj, assay = "SCT", npcs = 50)
    obj <- RunUMAP(obj, dims = 1:50)
