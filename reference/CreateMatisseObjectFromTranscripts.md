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

  A matrix or sparse matrix (**transcripts × cells**) of
  transcript-level counts. Row names must be transcript IDs matching
  those in the IOE files. Column names must be cell barcodes.

- ioe_files:

  Character vector of paths to SUPPA2 `.ioe` files. Multiple files (one
  per event type: SE, SS, MX, RI, FL) can be supplied.

- min_coverage:

  Integer. Minimum total transcript counts per cell per event required
  to report a PSI value. Default: `5`.

- verbose:

  Logical. Print progress messages. Default: `TRUE`.

## Value

A
[`MatisseObject`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
with `"transcript"` and `"psi"` assays populated inside the embedded
Seurat object.

## Details

Transcript counts are stored as a `Assay5` named `"transcript"` inside
the Seurat object (ready for
[`SCTransformTranscripts`](https://avisrilab.github.io/Matisse/reference/SCTransformTranscripts.md)).
PSI values, inclusion counts, and exclusion counts are stored inside a
`Assay5` named `"psi"` with layers `"data"`, `"counts"`, and
`"exclusion"`, respectively.

For each splice event the PSI per cell is: \$\$PSI\_{c,e} = \frac{\sum
\text{inclusion transcript counts}} {\sum \text{inclusion transcript
counts} + \sum \text{exclusion transcript counts}}\$\$

Cells where the total count (inclusion + exclusion) falls below
`min_coverage` are set to `NA`.

## See also

[`CreateMatisseObject`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md),
[`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md),
[`SCTransformTranscripts`](https://avisrilab.github.io/Matisse/reference/SCTransformTranscripts.md)
