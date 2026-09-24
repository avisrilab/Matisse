# Cluster cells in a MatisseObject

Runs
[`FindClusters`](https://satijalab.org/seurat/reference/FindClusters.html)
on the embedded Seurat object and returns the updated `MatisseObject`.
Cluster assignments are stored in `MatisseMeta(obj)$seurat_clusters` and
are immediately available for
[`PlotUMAP`](https://avisrilab.org/Matisse/reference/PlotUMAP.md) and
[`PlotViolin`](https://avisrilab.org/Matisse/reference/PlotViolin.md).

## Usage

``` r
# S3 method for class 'MatisseObject'
FindClusters(
  object,
  graph.name = NULL,
  cluster.name = NULL,
  modularity.fxn = 1,
  initial.membership = NULL,
  node.sizes = NULL,
  resolution = 0.8,
  method = NULL,
  algorithm = 1,
  leiden_method = c("leidenbase", "igraph"),
  leiden_objective_function = c("modularity", "CPM"),
  n.start = 10,
  n.iter = 10,
  random.seed = 0,
  group.singletons = TRUE,
  temp.file.location = NULL,
  edge.file.name = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A `MatisseObject`.

- graph.name:

  Name of graph to use for the clustering algorithm

- cluster.name:

  Name of output clusters

- modularity.fxn:

  Modularity function (1 = standard; 2 = alternative).

- initial.membership:

  Passed to the \`initial_membership\` parameter of
  \`leidenbase::leiden_find_partition\`.

- node.sizes:

  Passed to the \`node_sizes\` parameter of
  \`leidenbase::leiden_find_partition\`.

- resolution:

  Value of the resolution parameter, use a value above (below) 1.0 if
  you want to obtain a larger (smaller) number of communities.

- method:

  DEPRECATED.

- algorithm:

  Algorithm for modularity optimization (1 = original Louvain algorithm;
  2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4
  = Leiden algorithm).

- leiden_method:

  Choose from the leidenbase ("leidenbase") or igraph ("igraph")
  packages for running leiden. Default is "leidenbase"

- leiden_objective_function:

  objective function to use if \`leiden_method = "igraph"\`. See
  [`cluster_leiden`](https://r.igraph.org/reference/cluster_leiden.html)
  for more information. Default is "modularity".

- n.start:

  Number of random starts.

- n.iter:

  Maximal number of iterations per random start.

- random.seed:

  Seed of the random number generator.

- group.singletons:

  Group singletons into nearest cluster. If FALSE, assign all singletons
  to a "singleton" group

- temp.file.location:

  Directory where intermediate files will be written. Specify the
  ABSOLUTE path.

- edge.file.name:

  Edge file to use as input for modularity optimizer jar.

- verbose:

  Print output

- ...:

  Additional arguments forwarded to
  [`FindClusters`](https://satijalab.org/seurat/reference/FindClusters.html)
  (e.g. `resolution`).

## Value

The updated `MatisseObject` with `seurat_clusters` added to cell
metadata.

## See also

[`FindNeighbors.MatisseObject`](https://avisrilab.org/Matisse/reference/FindNeighbors.MatisseObject.md),
[`PlotUMAP`](https://avisrilab.org/Matisse/reference/PlotUMAP.md)
