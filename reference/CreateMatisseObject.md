# Create a MatisseObject

The primary constructor for
[`MatisseObject`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md).
Combines a `Seurat` object with optional isoform-resolved data layers.

## Usage

``` r
CreateMatisseObject(
  seurat,
  junction_counts = NULL,
  event_data = NULL,
  junction_data = NULL,
  verbose = TRUE
)
```

## Arguments

- seurat:

  A `Seurat` object. Required; provides cell barcodes, gene expression,
  and cell-level metadata.

- junction_counts:

  A sparse matrix (dgCMatrix, cells x junctions) of raw per-junction
  read counts. Row names must match `colnames(seurat)`. Default: `NULL`.

- event_data:

  A `data.frame` defining splice events. Required columns: `event_id`,
  `gene_id`, `chr`, `strand`, `event_type`, `inclusion_junctions`
  (semicolon-separated junction IDs), `exclusion_junctions`
  (semicolon-separated junction IDs). Default: `NULL` (empty table).

- junction_data:

  A `data.frame` of junction annotations. Required columns:
  `junction_id`, `chr`, `start`, `end`, `strand`, `gene_id`. Default:
  `NULL` (empty table).

- verbose:

  Logical. Print construction progress. Default: `TRUE`.

## Value

A
[`MatisseObject`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md).

## Examples

``` r
if (FALSE) { # \dontrun{
library(Seurat)
counts <- matrix(rpois(200, 5), nrow = 20,
                 dimnames = list(paste0("Gene", 1:20),
                                 paste0("Cell", 1:10)))
seu <- CreateSeuratObject(counts)
obj <- CreateMatisseObject(seurat = seu)
obj
} # }
```
