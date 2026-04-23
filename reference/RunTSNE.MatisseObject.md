# Run t-SNE on a MatisseObject

Runs [`RunTSNE`](https://satijalab.org/seurat/reference/RunTSNE.html) on
the embedded Seurat object and returns the updated `MatisseObject`. An
alternative to
[`RunUMAP.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunUMAP.MatisseObject.md)
for 2-D cell embedding.

## Usage

``` r
# S3 method for class 'MatisseObject'
RunTSNE(object, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments forwarded to
  [`RunTSNE`](https://satijalab.org/seurat/reference/RunTSNE.html) (e.g.
  `dims`, `perplexity`).

## Value

The updated `MatisseObject` with a `"tsne"` reduction.

## See also

[`RunUMAP.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunUMAP.MatisseObject.md),
[`RunPCA.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunPCA.MatisseObject.md)
