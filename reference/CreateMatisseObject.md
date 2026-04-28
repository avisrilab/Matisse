# Create a MatisseObject

The single constructor for
[`MatisseObject`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md).
Combines a `Seurat` object with isoform-resolved splicing data. The
operating mode is detected automatically from the inputs you supply:

## Usage

``` r
CreateMatisseObject(
  seurat,
  junction_counts = NULL,
  transcript_counts = NULL,
  ioe_files = NULL,
  event_data = NULL,
  junction_data = NULL,
  verbose = TRUE
)
```

## Arguments

- seurat:

  A `Seurat` object. Required.

- junction_counts:

  A sparse matrix (dgCMatrix, cells x junctions) of raw per-junction
  read counts. Row names must match `colnames(seurat)`. Triggers
  junction mode. Default: `NULL`.

- transcript_counts:

  A matrix or sparse matrix (transcripts x cells) of raw
  transcript-level counts. Stored as `Assay5("isoform")` in the Seurat
  object. Column names must overlap with `colnames(seurat)`. Triggers
  transcript mode. Default: `NULL`.

- ioe_files:

  Character vector of paths to SUPPA2 `.ioe` files. When supplied
  together with `transcript_counts`, the parsed event annotation is
  stored in the object and used by
  [`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md).
  Default: `NULL`.

- event_data:

  A `data.frame` defining splice events (junction mode only). Required
  columns: `event_id`, `gene_id`, `chr`, `strand`, `event_type`,
  `inclusion_junctions`, `exclusion_junctions`. Default: `NULL`.

- junction_data:

  A `data.frame` of junction annotations. Required columns:
  `junction_id`, `chr`, `start`, `end`, `strand`, `gene_id`. Default:
  `NULL`.

- verbose:

  Logical. Print construction progress. Default: `TRUE`.

## Value

A
[`MatisseObject`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md).

## Details

- **Junction mode** (short-read): pass `junction_counts`. Junction
  counts are stored as `Assay5("isoform")` inside the Seurat object.
  Call
  [`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
  afterwards to compute PSI values.

- **Transcript mode** (long-read): pass `transcript_counts` and
  optionally `ioe_files`. Transcript counts are stored as
  `Assay5("isoform")`. Call
  [`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
  (with `ioe_files`) to compute PSI values.

## Examples

``` r
if (FALSE) { # \dontrun{
library(Seurat)
counts <- matrix(rpois(200, 5), nrow = 20,
                 dimnames = list(paste0("Gene", 1:20),
                                 paste0("Cell", 1:10)))
seu <- CreateSeuratObject(counts)

# Junction mode
jxn <- make_junction_counts()
obj <- CreateMatisseObject(seurat = seu, junction_counts = jxn)
obj <- CalculatePSI(obj)

# Transcript mode
tx  <- make_transcript_counts()
obj <- CreateMatisseObject(seurat = seu, transcript_counts = tx,
                           ioe_files = "path/to/events.ioe")
obj <- CalculatePSI(obj)
} # }
```
