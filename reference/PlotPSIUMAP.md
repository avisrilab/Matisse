# UMAP plot colored by PSI of a specific splice event

Overlays the PSI value of a single splice event on the UMAP embedding
stored in the embedded Seurat object.

## Usage

``` r
PlotPSIUMAP(object, event_id, ...)

# S4 method for class 'MatisseObject'
PlotPSIUMAP(
  object,
  event_id,
  reduction = "umap",
  dims = c(1L, 2L),
  pt_size = 0.5,
  na_colour = "grey80",
  title = NULL
)
```

## Arguments

- object:

  A `MatisseObject` with a non-`NULL` `psi` slot. The Seurat object must
  contain a `umap` reduction.

- event_id:

  Character. Column name in the PSI matrix.

- ...:

  Additional arguments (see `PlotPSIUMAP`).

- reduction:

  Character. Name of the dimensionality reduction to use. Default:
  `"umap"`.

- dims:

  Integer vector of length 2 selecting which dimensions to plot.
  Default: `c(1, 2)`.

- pt_size:

  Numeric. Point size. Default: `0.5`.

- na_colour:

  Character. Colour for `NA` PSI cells. Default: `"grey80"`.

- title:

  Character. Plot title. Defaults to the event ID.

## Value

A `ggplot` object.

A `ggplot` object.
