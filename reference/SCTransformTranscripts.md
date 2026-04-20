# SCTransform normalization for the transcript assay

Runs `SCTransform` on the `"transcript"` assay of the embedded Seurat
object and follows it with `RunPCA` using a larger number of principal
components, giving better isoform-level cluster resolution.

Runs
[`SCTransform`](https://satijalab.org/seurat/reference/SCTransform.html)
on the `"transcript"` `Assay5` of the embedded Seurat object, followed
by [`RunPCA`](https://satijalab.org/seurat/reference/RunPCA.html) with
`n_pca_dims` components.

## Usage

``` r
SCTransformTranscripts(object, ...)

# S4 method for class 'MatisseObject'
SCTransformTranscripts(
  object,
  n_pca_dims = 50L,
  vars_to_regress = NULL,
  assay = "transcript",
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A
  [`MatisseObject`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
  with a `"transcript"` assay. Use
  [`CreateMatisseObjectFromTranscripts`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObjectFromTranscripts.md)
  or `CreateMatisseObject(transcript_counts = ...)` to create one.

- ...:

  Additional arguments forwarded to
  [`SCTransform`](https://satijalab.org/seurat/reference/SCTransform.html).

- n_pca_dims:

  Integer. Number of principal components to compute. Default: `50`.

- vars_to_regress:

  Character vector of variables to regress out during SCTransform (e.g.
  `"percent.mt"`). Default: `NULL`.

- assay:

  Character. Name of the assay to run SCTransform on. Default:
  `"transcript"`.

- verbose:

  Logical. Default: `TRUE`.

## Value

The updated `MatisseObject`.

The updated `MatisseObject` with:

- An `"SCT"` assay (created by SCTransform) in the Seurat object.

- A `"pca"` reduction computed on the SCT assay.

## Details

SCTransform is better suited than standard log-normalization for
transcript-level count data: it corrects for sequencing depth while
preserving biological variance. Using a higher number of PCA dimensions
(default 50) has been found to resolve more isoform-level clusters than
the standard 20.

After this step you can call
[`FindNeighbors`](https://satijalab.org/seurat/reference/FindNeighbors.html),
[`FindClusters`](https://satijalab.org/seurat/reference/FindClusters.html),
and [`RunUMAP`](https://satijalab.org/seurat/reference/RunUMAP.html)
directly on the `MatisseObject` via the built-in method forwarding (see
[`MatisseObject-class`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)).

## See also

[`CreateMatisseObjectFromTranscripts`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObjectFromTranscripts.md),
[`SCTransform`](https://satijalab.org/seurat/reference/SCTransform.html),
[`RunPCA`](https://satijalab.org/seurat/reference/RunPCA.html)
