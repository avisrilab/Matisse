# Scatter plot of a PSI dimensionality reduction

Plots cells in a 2-D PSI embedding coloured by a metadata variable.
Equivalent to
[`Seurat::DimPlot()`](https://satijalab.org/seurat/reference/DimPlot.html)
but operates on the PSI reductions stored in the embedded Seurat
(`"psi_umap"`, `"psi_pca"`, etc.).

## Usage

``` r
PlotPSIDimPlot(object, ...)

# S4 method for class 'MatisseObject'
PlotPSIDimPlot(
  object,
  reduction = "psi_umap",
  group_by = "psi_clusters",
  dims = c(1L, 2L),
  pt_size = 0.5,
  label = FALSE,
  title = NULL
)
```

## Arguments

- object:

  A `MatisseObject` whose embedded Seurat contains the reduction
  specified by `reduction`.

- ...:

  Additional arguments (see `PlotPSIDimPlot`).

- reduction:

  Character. Name of the reduction to plot. Default: `"psi_umap"`.

- group_by:

  Character. Column in `Seurat::meta.data` or
  [`MatisseMeta()`](https://avisrilab.github.io/Matisse/reference/MatisseMeta.md)
  used to colour cells. Default: `"psi_clusters"`.

- dims:

  Integer vector of length 2. Which dimensions to plot. Default:
  `c(1, 2)`.

- pt_size:

  Numeric. Point size. Default: `0.5`.

- label:

  Logical. Overlay cluster label at group centroid. Default: `FALSE`.

- title:

  Character. Plot title. Defaults to auto-generated.

## Value

A `ggplot` object.

A `ggplot` object.

## See also

[`RunPSIUMAP`](https://avisrilab.github.io/Matisse/reference/RunPSIUMAP.md),
[`FindPSIClusters`](https://avisrilab.github.io/Matisse/reference/FindPSIClusters.md)

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- RunPSIPCA(obj) |> RunPSIUMAP() |> FindPSIClusters()
PlotPSIDimPlot(obj, label = TRUE)
} # }
```
