# Run UMAP on a MatisseObject

Runs [`RunUMAP`](https://satijalab.org/seurat/reference/RunUMAP.html) on
the embedded Seurat object and returns the updated `MatisseObject`. Call
after
[`RunPCA.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunPCA.MatisseObject.md)
(or
[`RunSVD.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunSVD.MatisseObject.md)
for ATAC data). The resulting embedding is accessible via
`GetSeurat(obj)` and used by
[`PlotUMAP`](https://avisrilab.github.io/Matisse/reference/PlotUMAP.md).

## Usage

``` r
# S3 method for class 'MatisseObject'
RunUMAP(
  object,
  dims = NULL,
  reduction = "pca",
  features = NULL,
  graph = NULL,
  assay = DefaultAssay(object = object),
  nn.name = NULL,
  slot = "data",
  umap.method = "uwot",
  reduction.model = NULL,
  return.model = FALSE,
  n.neighbors = 30L,
  n.components = 2L,
  metric = "cosine",
  n.epochs = NULL,
  learning.rate = 1,
  min.dist = 0.3,
  spread = 1,
  set.op.mix.ratio = 1,
  local.connectivity = 1L,
  repulsion.strength = 1,
  negative.sample.rate = 5L,
  a = NULL,
  b = NULL,
  uwot.sgd = FALSE,
  uwot.approx_pow = FALSE,
  seed.use = 42L,
  metric.kwds = NULL,
  angular.rp.forest = FALSE,
  densmap = FALSE,
  dens.lambda = 2,
  dens.frac = 0.3,
  dens.var.shift = 0.1,
  verbose = TRUE,
  reduction.name = "umap",
  reduction.key = NULL,
  ...
)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments forwarded to
  [`RunUMAP`](https://satijalab.org/seurat/reference/RunUMAP.html) (e.g.
  `dims`, `reduction`).

## Value

The updated `MatisseObject` with a `"umap"` reduction.

## See also

[`RunPCA.MatisseObject`](https://avisrilab.github.io/Matisse/reference/RunPCA.MatisseObject.md),
[`PlotUMAP`](https://avisrilab.github.io/Matisse/reference/PlotUMAP.md)
