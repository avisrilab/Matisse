# Create a MatisseObject from transcript-level counts

Constructs a `MatisseObject` directly from a transcript-by-cell count
matrix (e.g. from Bagpiper or other transcript-aware quantifiers) by
mapping transcripts to splice events using SUPPA2 IOE annotation files.

## Usage

``` r
CreateMatisseObjectFromTranscripts(
  seurat,
  transcript_counts,
  ioe_files,
  min_coverage = 5L,
  verbose = TRUE
)
```

## Arguments

- seurat:

  A `Seurat` object. Provides cell barcodes and gene-level metadata.

- transcript_counts:

  A matrix or sparse matrix (**transcripts x cells**) of
  transcript-level counts. Row names must be transcript IDs matching
  those in the IOE files (e.g. GENCODE/Ensembl IDs). Column names must
  be cell barcodes.

- ioe_files:

  Character vector of paths to SUPPA2 `.ioe` files. Multiple files (one
  per event type: SE, SS, MX, RI, FL) can be supplied and will be
  combined. Each file must have the standard four-column SUPPA2 IOE
  format: `seqname`, `gene_id;event_id`, `inclusion_transcripts`,
  `total_transcripts`.

- min_coverage:

  Integer. Minimum total transcript counts (inclusion + exclusion) per
  cell per event required to report a PSI value. Cells below this
  threshold receive `NA`. Default: `5`.

- verbose:

  Logical. Print progress messages. Default: `TRUE`.

## Value

A
[`MatisseObject`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
with `psi`, `inclusion_counts`, and `exclusion_counts` slots populated.
The `junction_counts` slot is `NULL` (not applicable for
transcript-based input). The `event_data` slot uses transcript IDs in
the `inclusion_junctions` and `exclusion_junctions` columns.

## Details

For each splice event the PSI per cell is: \$\$PSI\_{c,e} = \frac{\sum
\text{inclusion transcript counts}} {\sum \text{inclusion transcript
counts} + \sum \text{exclusion transcript counts}}\$\$

Cells where the total count (inclusion + exclusion) falls below
`min_coverage` are set to `NA`.

## See also

[`CreateMatisseObject`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md),
[`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming you have:
#   - a Seurat object `seu`
#   - a transcripts-x-cells matrix `tx_mat` with GENCODE transcript IDs
#   - SUPPA2 IOE files for SE and RI events
obj <- CreateMatisseObjectFromTranscripts(
  seurat            = seu,
  transcript_counts = tx_mat,
  ioe_files         = c("events_SE.ioe", "events_RI.ioe"),
  min_coverage      = 5L
)
} # }
```
