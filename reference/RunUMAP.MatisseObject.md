# Run UMAP on a MatisseObject

Runs [`RunUMAP`](https://satijalab.org/seurat/reference/RunUMAP.html) on
the embedded Seurat object and returns the updated `MatisseObject`. Call
after
[`RunPCA.MatisseObject`](https://avisrilab.org/Matisse/reference/RunPCA.MatisseObject.md)
(or
[`RunSVD.MatisseObject`](https://avisrilab.org/Matisse/reference/RunSVD.MatisseObject.md)
for ATAC data). The resulting embedding is accessible via
`GetSeurat(obj)` and used by
[`PlotUMAP`](https://avisrilab.org/Matisse/reference/PlotUMAP.md).

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
  uwot.init = "spectral",
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

- dims:

  Which dimensions to use as input features, used only if `features` is
  NULL

- reduction:

  Which dimensional reduction (PCA or ICA) to use for the UMAP input.
  Default is PCA

- features:

  If set, run UMAP on this subset of features (instead of running on a
  set of reduced dimensions). Not set (NULL) by default; `dims` must be
  NULL to run on features

- graph:

  Name of graph on which to run UMAP

- assay:

  Assay to pull data for when using `features`, or assay used to
  construct Graph if running UMAP on a Graph

- nn.name:

  Name of knn output on which to run UMAP

- slot:

  The slot used to pull data for when using `features`. data slot is by
  default.

- umap.method:

  UMAP implementation to run. Can be

  `uwot`:

  :   Runs umap via the uwot R package
      [`umap`](https://jlmelville.github.io/uwot/reference/umap.html)

  `uwot2`:

  :   Runs umap2 via the uwot R package
      [`umap2`](https://jlmelville.github.io/uwot/reference/umap2.html)

  `umap-learn`:

  :   Run the Seurat wrapper of the python umap-learn package

- reduction.model:

  `DimReduc` object that contains the umap model

- return.model:

  whether UMAP will return the uwot model

- n.neighbors:

  This determines the number of neighboring points used in local
  approximations of manifold structure. Larger values will result in
  more global structure being preserved at the loss of detailed local
  structure. In general this parameter should often be in the range 5 to
  50.

- n.components:

  The dimension of the space to embed into.

- metric:

  metric: This determines the choice of metric used to measure distance
  in the input space. A wide variety of metrics are already coded, and a
  user defined function can be passed as long as it has been JITd by
  numba.

- n.epochs:

  he number of training epochs to be used in optimizing the low
  dimensional embedding. Larger values result in more accurate
  embeddings. If NULL is specified, a value will be selected based on
  the size of the input dataset (200 for large datasets, 500 for small).

- learning.rate:

  The initial learning rate for the embedding optimization.

- min.dist:

  This controls how tightly the embedding is allowed compress points
  together. Larger values ensure embedded points are more evenly
  distributed, while smaller values allow the algorithm to optimize more
  accurately with regard to local structure. Sensible values are in the
  range 0.001 to 0.5.

- spread:

  The effective scale of embedded points. In combination with min.dist
  this determines how clustered/clumped the embedded points are.

- set.op.mix.ratio:

  Interpolate between (fuzzy) union and intersection as the set
  operation used to combine local fuzzy simplicial sets to obtain a
  global fuzzy simplicial sets. Both fuzzy set operations use the
  product t-norm. The value of this parameter should be between 0.0 and
  1.0; a value of 1.0 will use a pure fuzzy union, while 0.0 will use a
  pure fuzzy intersection.

- local.connectivity:

  The local connectivity required - i.e. the number of nearest neighbors
  that should be assumed to be connected at a local level. The higher
  this value the more connected the manifold becomes locally. In
  practice this should be not more than the local intrinsic dimension of
  the manifold.

- repulsion.strength:

  Weighting applied to negative samples in low dimensional embedding
  optimization. Values higher than one will result in greater weight
  being given to negative samples.

- negative.sample.rate:

  The number of negative samples to select per positive sample in the
  optimization process. Increasing this value will result in greater
  repulsive force being applied, greater optimization cost, but slightly
  more accuracy.

- a:

  More specific parameters controlling the embedding. If NULL, these
  values are set automatically as determined by min. dist and spread.
  Parameter of differentiable approximation of right adjoint functor.

- b:

  More specific parameters controlling the embedding. If NULL, these
  values are set automatically as determined by min. dist and spread.
  Parameter of differentiable approximation of right adjoint functor.

- uwot.sgd:

  Set `uwot::umap(fast_sgd = TRUE)`; see
  [`umap`](https://jlmelville.github.io/uwot/reference/umap.html) for
  more details

- uwot.approx_pow:

  Set `uwot::umap(approx_pow = TRUE)`. Default is `FALSE`. See
  [`umap`](https://jlmelville.github.io/uwot/reference/umap.html) for
  more details.

- uwot.init:

  Set the initialization method to use for
  [`umap`](https://jlmelville.github.io/uwot/reference/umap.html) or
  [`umap2`](https://jlmelville.github.io/uwot/reference/umap2.html),
  which is passed to the `init` parameter of these functions. See these
  functions for available options. Default is `"spectral"`.

- seed.use:

  Set a random seed. By default, sets the seed to 42. Setting NULL will
  not set a seed

- metric.kwds:

  A dictionary of arguments to pass on to the metric, such as the p
  value for Minkowski distance. If NULL then no arguments are passed on.

- angular.rp.forest:

  Whether to use an angular random projection forest to initialize the
  approximate nearest neighbor search. This can be faster, but is mostly
  on useful for metric that use an angular style distance such as
  cosine, correlation etc. In the case of those metrics angular forests
  will be chosen automatically.

- densmap:

  Whether to use the density-augmented objective of densMAP. Turning on
  this option generates an embedding where the local densities are
  encouraged to be correlated with those in the original space.
  Parameters below with the prefix ‘dens’ further control the behavior
  of this extension. Default is FALSE. Only compatible with 'umap-learn'
  method and version of umap-learn \>= 0.5.0

- dens.lambda:

  Specific parameter which controls the regularization weight of the
  density correlation term in densMAP. Higher values prioritize density
  preservation over the UMAP objective, and vice versa for values closer
  to zero. Setting this parameter to zero is equivalent to running the
  original UMAP algorithm. Default value is 2.

- dens.frac:

  Specific parameter which controls the fraction of epochs (between 0
  and 1) where the density-augmented objective is used in densMAP. The
  first (1 - dens_frac) fraction of epochs optimize the original UMAP
  objective before introducing the density correlation term. Default is
  0.3.

- dens.var.shift:

  Specific parameter which specifies a small constant added to the
  variance of local radii in the embedding when calculating the density
  correlation objective to prevent numerical instability from dividing
  by a small number. Default is 0.1.

- verbose:

  Controls verbosity

- reduction.name:

  Name to store dimensional reduction under in the Seurat object

- reduction.key:

  dimensional reduction key, specifies the string before the number for
  the dimension names. UMAP by default

- ...:

  Additional arguments forwarded to
  [`RunUMAP`](https://satijalab.org/seurat/reference/RunUMAP.html) (e.g.
  `dims`, `reduction`).

## Value

The updated `MatisseObject` with a `"umap"` reduction.

## See also

[`RunPCA.MatisseObject`](https://avisrilab.org/Matisse/reference/RunPCA.MatisseObject.md),
[`PlotUMAP`](https://avisrilab.org/Matisse/reference/PlotUMAP.md)
