# Find clusters in PSI space

Builds a shared-nearest-neighbour (SNN) graph on the PSI PCA embedding
via
[`Seurat::FindNeighbors()`](https://satijalab.org/seurat/reference/FindNeighbors.html),
then identifies communities using the Louvain algorithm via
[`Seurat::FindClusters()`](https://satijalab.org/seurat/reference/FindClusters.html).
Cluster labels are stored as `psi_clusters` in `object@seurat@meta.data`
and are accessible as `object\$psi_clusters`.

## Usage

``` r
FindPSIClusters(object, ...)

# S4 method for class 'MatisseObject'
FindPSIClusters(
  object,
  resolution = 0.8,
  dims = 1:20,
  reduction = "psi_pca",
  verbose = TRUE
)
```

## Arguments

- object:

  A `MatisseObject` whose embedded Seurat contains the reduction
  specified by `reduction`.

- ...:

  Additional arguments (see `FindPSIClusters`).

- resolution:

  Numeric. Louvain resolution parameter (higher value → more clusters).
  Default: `0.8`.

- dims:

  Integer vector. PCA dimensions to use for the SNN graph. Default:
  `1:20`.

- reduction:

  Character. Name of the input reduction. Default: `"psi_pca"`.

- verbose:

  Logical. Default: `TRUE`.

## Value

The updated `MatisseObject` with `psi_clusters` metadata.

The input `MatisseObject` with `psi_clusters` added to
`object@seurat@meta.data`.

## See also

[`RunPSIPCA`](https://avisrilab.github.io/Matisse/reference/RunPSIPCA.md),
[`PlotPSIDimPlot`](https://avisrilab.github.io/Matisse/reference/PlotPSIDimPlot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- RunPSIPCA(obj, n_pcs = 30)
obj <- FindPSIClusters(obj, resolution = 0.5, dims = 1:20)
table(obj$psi_clusters)
} # }
```
