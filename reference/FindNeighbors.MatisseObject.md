# Compute a shared nearest-neighbour graph for a MatisseObject

Runs
[`FindNeighbors`](https://satijalab.org/seurat/reference/FindNeighbors.html)
on the embedded Seurat object and returns the updated `MatisseObject`.
Typically called after
[`RunPCA.MatisseObject`](https://avisrilab.org/Matisse/reference/RunPCA.MatisseObject.md)
and before
[`FindClusters.MatisseObject`](https://avisrilab.org/Matisse/reference/FindClusters.MatisseObject.md).

## Usage

``` r
# S3 method for class 'MatisseObject'
FindNeighbors(
  object,
  reduction = "pca",
  dims = 1:10,
  assay = NULL,
  features = NULL,
  k.param = 20,
  return.neighbor = FALSE,
  compute.SNN = !return.neighbor,
  prune.SNN = 1/15,
  nn.method = "annoy",
  n.trees = 50,
  annoy.metric = "euclidean",
  nn.eps = 0,
  verbose = TRUE,
  do.plot = FALSE,
  graph.name = NULL,
  l2.norm = FALSE,
  cache.index = FALSE,
  ...
)
```

## Arguments

- object:

  A `MatisseObject`.

- reduction:

  Reduction to use as input for building the (S)NN

- dims:

  Dimensions of reduction to use as input

- assay:

  Assay to use in construction of (S)NN; used only when `dims` is `NULL`

- features:

  Features to use as input for building the (S)NN; used only when `dims`
  is `NULL`

- k.param:

  Defines k for the k-nearest neighbor algorithm

- return.neighbor:

  Return result as
  [`Neighbor`](https://satijalab.github.io/seurat-object/reference/Neighbor-class.html)
  object. Not used with distance matrix input.

- compute.SNN:

  also compute the shared nearest neighbor graph

- prune.SNN:

  Sets the cutoff for acceptable Jaccard index when computing the
  neighborhood overlap for the SNN construction. Any edges with values
  less than or equal to this will be set to 0 and removed from the SNN
  graph. Essentially sets the stringency of pruning (0 — no pruning, 1 —
  prune everything).

- nn.method:

  Method for nearest neighbor finding. Options include: rann, annoy

- n.trees:

  More trees gives higher precision when using annoy approximate nearest
  neighbor search

- annoy.metric:

  Distance metric for annoy. Options include: euclidean, cosine,
  manhattan, and hamming

- nn.eps:

  Error bound when performing nearest neighbor search using RANN;
  default of 0.0 implies exact nearest neighbor search

- verbose:

  Whether or not to print output to the console

- do.plot:

  Plot SNN graph on tSNE coordinates

- graph.name:

  Optional naming parameter for stored (S)NN graph (or Neighbor object,
  if return.neighbor = TRUE). Default is assay.name\_(s)nn. To store
  both the neighbor graph and the shared nearest neighbor (SNN) graph,
  you must supply a vector containing two names to the `graph.name`
  parameter. The first element in the vector will be used to store the
  nearest neighbor (NN) graph, and the second element used to store the
  SNN graph. If only one name is supplied, only the NN graph is stored.

- l2.norm:

  Take L2Norm of the data

- cache.index:

  Include cached index in returned Neighbor object (only relevant if
  return.neighbor = TRUE)

- ...:

  Additional arguments forwarded to
  [`FindNeighbors`](https://satijalab.org/seurat/reference/FindNeighbors.html)
  (e.g. `dims`, `reduction`).

## Value

The updated `MatisseObject` with neighbour graph stored inside the
embedded Seurat object.

## See also

[`FindClusters.MatisseObject`](https://avisrilab.org/Matisse/reference/FindClusters.MatisseObject.md),
[`RunPCA.MatisseObject`](https://avisrilab.org/Matisse/reference/RunPCA.MatisseObject.md)
