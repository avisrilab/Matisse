# Run UMAP on PSI principal components

Computes a 2-D UMAP embedding from a PSI PCA reduction (or any other
named reduction) using the uwot package. The result is stored as
`"psi_umap"` in the embedded Seurat object.

## Usage

``` r
RunPSIUMAP(object, ...)

# S4 method for class 'MatisseObject'
RunPSIUMAP(
  object,
  dims = 1:30,
  reduction = "psi_pca",
  n_neighbors = 30L,
  min_dist = 0.3,
  seed_use = 42L,
  verbose = TRUE
)
```

## Arguments

- object:

  A `MatisseObject` whose embedded Seurat contains the reduction
  specified by `reduction` (typically after calling
  [`RunPSIPCA`](https://avisrilab.github.io/Matisse/reference/RunPSIPCA.md)).

- ...:

  Additional arguments (see `RunPSIUMAP`).

- dims:

  Integer vector. Which PCA dimensions to use as UMAP input. Default:
  `1:30`.

- reduction:

  Character. Name of the input reduction. Default: `"psi_pca"`.

- n_neighbors:

  Integer. UMAP `n_neighbors` parameter. Default: `30`.

- min_dist:

  Numeric. UMAP `min_dist` parameter. Default: `0.3`.

- seed_use:

  Integer. Random seed for reproducibility. Default: `42`.

- verbose:

  Logical. Default: `TRUE`.

## Value

The updated `MatisseObject` with `"psi_umap"` reduction.

The input `MatisseObject` with `"psi_umap"` added to the embedded Seurat
reductions.

## See also

[`RunPSIPCA`](https://avisrilab.github.io/Matisse/reference/RunPSIPCA.md),
[`FindPSIClusters`](https://avisrilab.github.io/Matisse/reference/FindPSIClusters.md),
[`PlotPSIDimPlot`](https://avisrilab.github.io/Matisse/reference/PlotPSIDimPlot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- RunPSIPCA(obj, n_pcs = 30)
obj <- RunPSIUMAP(obj, dims = 1:20)
} # }
```
