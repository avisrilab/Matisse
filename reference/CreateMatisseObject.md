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
  min_coverage = 5L,
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
  transcript-level counts. Stored as `Assay5("transcript")` in the
  Seurat object. Column names must overlap with `colnames(seurat)`.
  Default: `NULL`.

- ioe_files:

  Character vector of paths to SUPPA2 `.ioe` files. When supplied
  together with `transcript_counts`, triggers event mode and PSI is
  computed at construction. Default: `NULL`.

- event_data:

  A `data.frame` defining splice events (junction mode only). Required
  columns: `event_id`, `gene_id`, `chr`, `strand`, `event_type`,
  `inclusion_junctions`, `exclusion_junctions`. Default: `NULL`.

- junction_data:

  A `data.frame` of junction annotations. Required columns:
  `junction_id`, `chr`, `start`, `end`, `strand`, `gene_id`. Default:
  `NULL`.

- min_coverage:

  Integer. Minimum total transcript counts per cell per event to report
  a PSI value (event mode only). Default: `5`.

- verbose:

  Logical. Print construction progress. Default: `TRUE`.

## Value

A
[`MatisseObject`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md).

## Details

- **Junction mode** (short-read): pass `junction_counts`. Junction
  counts are stored as `Assay5("junction")` inside the Seurat object.
  Call
  [`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
  afterwards to compute PSI values.

- **Event mode** (long-read): pass `transcript_counts` and `ioe_files`.
  Transcript counts are stored as `Assay5("transcript")` and PSI is
  computed immediately from the supplied SUPPA2 `.ioe` event
  definitions, stored as `Assay5("psi")`.

You can also pass `transcript_counts` alone (without `ioe_files`) to
store the transcript assay without computing PSI – for example, when you
want to run
[`SCTransform`](https://satijalab.org/seurat/reference/SCTransform.html)
on transcript-level counts and compute PSI separately.

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

# Event mode
tx  <- make_transcript_counts()
obj <- CreateMatisseObject(seurat = seu, transcript_counts = tx,
                           ioe_files = "path/to/events.ioe")
} # }
```
