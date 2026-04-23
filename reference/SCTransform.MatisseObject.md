# SCTransform normalisation for MatisseObjects

Runs
[`SCTransform`](https://satijalab.org/seurat/reference/SCTransform.html)
with mode-aware defaults. In **event mode** (long-read), normalises the
`"transcript"` assay by default, so that transcript-level abundances are
variance-stabilised before dimensionality reduction. In **junction
mode** (short-read), normalises the active default assay (typically
`"RNA"`).

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

  Character vector. Covariates to regress out (e.g. `"percent.mt"`).
  Default: `NULL`.

- verbose:

  Logical. Print progress messages. Default: `TRUE`.

- ...:

  Additional arguments forwarded to
  [`SCTransform`](https://satijalab.org/seurat/reference/SCTransform.html).

## Value

The updated `MatisseObject` with a new `"SCT"` assay.

## Details

Override the target assay with the `assay` argument. PCA is **not** run
automatically — call
[`RunPCA.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunPCA.MatisseObject.md)
after normalisation:

    obj <- SCTransform(obj)                        # normalise
    obj <- RunPCA(obj, assay = "SCT", npcs = 50)  # reduce
    obj <- RunUMAP(obj, dims = 1:50)              # embed

## See also

[`RunPCA.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunPCA.MatisseObject.md),
[`NormalizeData.MatisseObject`](https://avisrilab.github.io/Matisse/reference/NormalizeData.MatisseObject.md)
