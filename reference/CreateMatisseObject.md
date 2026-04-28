# Create a MatisseObject

The single constructor for
[`MatisseObject`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md).
Combines a `Seurat` object with isoform-resolved splicing data and
*also* computes PSI in one step (unless `defer_psi = TRUE`). The
operating mode is detected automatically from the inputs you supply:

## Usage

``` r
CreateMatisseObject(
  seurat,
  junction_counts = NULL,
  transcript_counts = NULL,
  events = NULL,
  min_coverage = 5L,
  defer_psi = FALSE,
  verbose = TRUE
)
```

## Arguments

- seurat:

  A `Seurat` object. Required.

- junction_counts:

  A sparse matrix (dgCMatrix, cells x junctions) of raw per-junction
  read counts. Row names must match `colnames(seurat)`. Junction IDs
  (column names) ideally encode genomic coordinates
  (`chr-start-end-strand`, `chr:start-end:strand`, or
  `chr_start_end_strand`); Matisse parses these for sashimi plots.
  Triggers junction mode. Default: `NULL`.

- transcript_counts:

  A matrix or sparse matrix (transcripts x cells) of raw
  transcript-level counts. Stored as `Assay5("isoform")` in the Seurat
  object. Column names must overlap with `colnames(seurat)`. Triggers
  transcript mode. Default: `NULL`.

- events:

  Splice-event annotation. Either:

  - a character vector of paths to SUPPA2 `.ioe` file(s) (parsed
    internally), or

  - a pre-built `data.frame` with columns `event_id`, `gene_id`, `chr`,
    `strand`, `event_type`, `inclusion_features`, `exclusion_features`.
    The `*_features` columns hold the IDs that support each role —
    junction IDs in junction mode, transcript IDs in transcript mode.

  Required.

- min_coverage:

  Integer. Minimum total reads per cell per event for the PSI
  calculation step. Default: `5L`. Forwarded to
  [`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
  when `defer_psi = FALSE`.

- defer_psi:

  Logical. Skip the PSI calculation step at construction. The returned
  object will have raw counts only (no `"psi"` assay). Call
  [`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
  later. Default: `FALSE`.

- verbose:

  Logical. Print construction progress. Default: `TRUE`.

## Value

A
[`MatisseObject`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md).

## Details

- **Junction mode** (short-read): pass `junction_counts` plus splice
  events via `events`. Junction counts are stored as `Assay5("isoform")`
  inside the Seurat object; PSI values land in `Assay5("psi")`.

- **Transcript mode** (long-read): pass `transcript_counts` plus splice
  events via `events` (typically a SUPPA2 IOE file). Transcript counts
  are stored as `Assay5("isoform")`; PSI values land in `Assay5("psi")`.

Matisse intentionally does **not** ship a heuristic event caller — bring
your own splice-event annotation (SUPPA2 IOE, rMATS table, or
equivalent).

## Examples

``` r
if (FALSE) { # \dontrun{
library(Seurat)
counts <- matrix(rpois(200, 5), nrow = 20,
                 dimnames = list(paste0("Gene", 1:20),
                                 paste0("Cell", 1:10)))
seu <- CreateSeuratObject(counts)

# Junction mode (short-read): events as a data.frame
jxn <- make_junction_counts()
obj <- CreateMatisseObject(seurat = seu, junction_counts = jxn,
                           events = my_event_df)

# Transcript mode (long-read): events as SUPPA2 IOE file paths
tx  <- make_transcript_counts()
obj <- CreateMatisseObject(seurat = seu, transcript_counts = tx,
                           events = "path/to/events.ioe")
} # }
```
