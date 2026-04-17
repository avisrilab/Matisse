# Run PCA on the PSI matrix

Imputes missing PSI values (uncovered entries) with per-event column
means, optionally filters low-variance events, then runs principal
component analysis via
[`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html). The result is
stored as a `DimReduc` named `"psi_pca"` inside the embedded Seurat
object and is accessible via
`SeuratObject::Embeddings(GetSeurat(object), "psi_pca")`.

## Usage

``` r
RunPSIPCA(object, ...)

# S4 method for class 'MatisseObject'
RunPSIPCA(
  object,
  n_pcs = 30L,
  events = NULL,
  cells = NULL,
  min_variance = NULL,
  center = TRUE,
  scale = FALSE,
  verbose = TRUE
)
```

## Arguments

- object:

  A `MatisseObject` with a non-`NULL` `psi` slot.

- ...:

  Additional arguments (see `RunPSIPCA`).

- n_pcs:

  Integer. Number of principal components to compute. Default: `30`.

- events:

  Character vector of event IDs to include. Default: all events.

- cells:

  Character vector of cell barcodes to include. Default: all cells.

- min_variance:

  Numeric. Retain only events whose variance (after NA imputation) is
  \>= this value. Default: `NULL` (no filtering).

- center:

  Logical. Center the PSI matrix before PCA. Default: `TRUE`.

- scale:

  Logical. Scale each event to unit variance. Default: `FALSE`.

- verbose:

  Logical. Print progress messages. Default: `TRUE`.

## Value

The updated `MatisseObject` with `"psi_pca"` reduction.

The input `MatisseObject` with `"psi_pca"` added to the embedded Seurat
reductions.

## See also

[`RunPSIUMAP`](https://avisrilab.github.io/Matisse/reference/RunPSIUMAP.md),
[`FindPSIClusters`](https://avisrilab.github.io/Matisse/reference/FindPSIClusters.md),
[`PlotPSIDimPlot`](https://avisrilab.github.io/Matisse/reference/PlotPSIDimPlot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- RunPSIPCA(obj, n_pcs = 30)
SeuratObject::Embeddings(GetSeurat(obj), "psi_pca")[1:5, 1:3]
} # }
```
