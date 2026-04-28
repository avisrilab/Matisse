# UMAP plot coloured by any feature

Overlays the value of a feature on the UMAP embedding stored in the
embedded Seurat object. Pass a PSI event ID, a junction ID, a gene name,
or a cell-metadata column. The colour scale adapts automatically:
diverging RdBu \[0,1\] for PSI; sequential viridis for counts and
expression.

Overlays the value of a feature on the UMAP embedding stored in the
embedded Seurat object. The colour scale adapts to the feature type:
diverging RdBu \[0,1\] for PSI events; sequential viridis for junction
counts and gene expression.

## Usage

``` r
PlotUMAP(
  object,
  feature,
  reduction = "umap",
  dims = c(1L, 2L),
  pt_size = 0.5,
  na_colour = "grey80",
  title = NULL,
  ...
)

# S4 method for class 'MatisseObject'
PlotUMAP(
  object,
  feature,
  reduction = "umap",
  dims = c(1L, 2L),
  pt_size = 0.5,
  na_colour = "grey80",
  title = NULL,
  ...
)
```

## Arguments

- object:

  A `MatisseObject` with a UMAP reduction.

- feature:

  Character. Feature to plot. May be a PSI event ID (e.g.
  `"SE:chr1:100-200:300-400:+"`), a junction or transcript ID, a gene
  name, or a cell-metadata column.

- reduction:

  Character. Name of the dimensionality reduction to use. Default:
  `"umap"`.

- dims:

  Integer vector of length 2 selecting which dimensions to plot.
  Default: `c(1, 2)`.

- pt_size:

  Numeric. Point size. Default: `0.5`.

- na_colour:

  Character. Colour for cells with no data. Default: `"grey80"`.

- title:

  Character. Plot title. Defaults to the feature name.

- ...:

  Additional arguments (see `PlotUMAP`).

## Value

A `ggplot` object.

A `ggplot` object.
