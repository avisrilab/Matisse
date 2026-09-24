# Find differentially expressed markers for a MatisseObject

Runs
[`FindMarkers`](https://satijalab.org/seurat/reference/FindMarkers.html)
on the embedded Seurat object. Unlike most dispatch methods, this
returns a `data.frame` of marker statistics rather than an updated
`MatisseObject`.

## Usage

``` r
# S3 method for class 'MatisseObject'
FindMarkers(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
  latent.vars = NULL,
  group.by = NULL,
  subset.ident = NULL,
  assay = NULL,
  reduction = NULL,
  ...
)
```

## Arguments

- object:

  A `MatisseObject`.

- ident.1:

  Identity class to define markers for; pass an object of class `phylo`
  or 'clustertree' to find markers for a node in a cluster tree; passing
  'clustertree' requires
  [`BuildClusterTree`](https://satijalab.org/seurat/reference/BuildClusterTree.html)
  to have been run

- ident.2:

  A second identity class for comparison; if `NULL`, use all other cells
  for comparison; if an object of class `phylo` or 'clustertree' is
  passed to `ident.1`, must pass a node to find markers for

- latent.vars:

  Variables to test, used only when `test.use` is one of 'LR',
  'negbinom', 'poisson', or 'MAST'

- group.by:

  Regroup cells into a different identity class prior to performing
  differential expression (see example); `"ident"` to use Idents

- subset.ident:

  Subset a particular identity class prior to regrouping. Only relevant
  if group.by is set (see example)

- assay:

  Assay to use in differential expression testing

- reduction:

  Reduction to use in differential expression testing - will test for DE
  on cell embeddings

- ...:

  Additional arguments forwarded to
  [`FindMarkers`](https://satijalab.org/seurat/reference/FindMarkers.html)
  (e.g. `ident.1`, `ident.2`, `group.by`, `features`).

## Value

A `data.frame` of marker genes with columns `p_val`, `avg_log2FC`,
`pct.1`, `pct.2`, `p_val_adj`.

## See also

[`FindClusters.MatisseObject`](https://avisrilab.org/Matisse/reference/FindClusters.MatisseObject.md)
