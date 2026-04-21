# SCTransform normalisation for MatisseObjects

Runs
[`SCTransform`](https://satijalab.org/seurat/reference/SCTransform.html)
with mode-aware defaults. In **event mode**, normalises the
`"transcript"` assay. In **junction mode**, normalises the active
default assay (usually `"RNA"`). Override with the `assay` argument.

## Usage

``` r
# S3 method for class 'MatisseObject'
SCTransform(
  object,
  assay = NULL,
  n_pca_dims = 50L,
  vars_to_regress = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A `MatisseObject`.

- assay:

  Character. Assay to normalise. Default: `"transcript"` in event mode;
  the active default assay in junction mode.

- n_pca_dims:

  Integer. PCA dimensions to compute. Default: `50`. Set to `0L` to skip
  PCA.

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

After normalisation, `RunPCA` is run automatically on the resulting
`"SCT"` assay. To skip PCA set `n_pca_dims = 0L`.
