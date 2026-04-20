# Create a MatisseObject

The primary constructor for
[`MatisseObject`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md).
Combines a `Seurat` object with optional isoform-resolved data layers.
When `transcript_counts` is supplied the matrix is stored as a `Assay5`
named `"transcript"` inside the embedded Seurat object, ready for
downstream SCTransform normalisation.

## Usage

``` r
CreateMatisseObject(
  seurat,
  junction_counts = NULL,
  transcript_counts = NULL,
  event_data = NULL,
  junction_data = NULL,
  verbose = TRUE
)
```

## Arguments

- seurat:

  A `Seurat` object. Required.

- junction_counts:

  A sparse matrix (dgCMatrix, cells × junctions) of raw per-junction
  read counts. Row names must match `colnames(seurat)`. Default: `NULL`.

- transcript_counts:

  A matrix or sparse matrix (transcripts × cells) of raw
  transcript-level counts. When supplied it is added to the Seurat
  object as a `Assay5` named `"transcript"`. Column names must overlap
  with `colnames(seurat)`. Default: `NULL`.

- event_data:

  A `data.frame` defining splice events. Required columns: `event_id`,
  `gene_id`, `chr`, `strand`, `event_type`, `inclusion_junctions`,
  `exclusion_junctions`. Default: `NULL`.

- junction_data:

  A `data.frame` of junction annotations. Required columns:
  `junction_id`, `chr`, `start`, `end`, `strand`, `gene_id`. Default:
  `NULL`.

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
