# Read a long-read transcript count matrix

Loads a per-cell transcript count matrix from a 10x-style MatrixMarket
triplet (`matrix.mtx`, `barcodes.tsv`, `features.tsv`, optionally `.gz`)
as produced by long-read / isoform quantifiers such as **Bagpiper**,
**FLAMES**, or **LIQA**, and returns it in the *transcripts x cells*
orientation with dimnames that
[`CreateMatisseObject`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md)'s
transcript mode expects.

## Usage

``` r
ReadTranscriptMatrix(
  mtx_dir,
  cells = NULL,
  strip_version = TRUE,
  orientation = c("auto", "cells_x_tx", "tx_x_cells"),
  verbose = TRUE
)
```

## Arguments

- mtx_dir:

  Directory holding `matrix.mtx`, `barcodes.tsv` and `features.tsv`
  (optionally `.gz`).

- cells:

  Optional character vector of cell barcodes to keep (e.g. filtered
  barcodes). Barcodes absent from the matrix are dropped with a warning.
  Default: `NULL` (keep all).

- strip_version:

  Logical. Strip a trailing `.<digits>` version suffix from transcript
  IDs; counts of IDs that collapse to the same stripped ID are summed.
  Default: `TRUE`.

- orientation:

  One of `"auto"` (infer from dimensions; default), `"cells_x_tx"`
  (source is cells x transcripts; transpose), or `"tx_x_cells"` (source
  is already transcripts x cells).

- verbose:

  Logical. Print progress. Default: `TRUE`.

## Value

A `dgCMatrix` of transcript counts, *transcripts x cells*, transcript
IDs as row names and barcodes as column names. Feed directly to
`CreateMatisseObject(transcript_counts = ...)`.

## Details

Quantifiers disagree on matrix orientation (Bagpiper writes *cells x
transcripts*; others write *transcripts x cells*), so
`orientation = "auto"` infers it by matching the two matrix dimensions
against the feature and barcode counts and transposes as needed.

**Transcript-version landmine.** SUPPA2 `.ioe` events and the quantifier
may disagree on whether transcript IDs carry a version suffix
(`ENST...\.3`). PSI aggregation (`.build_indicator_matrix`) matches IDs
only after a `_`-\>`-` sanitisation; it never strips versions, so a
mismatch silently zeroes PSI. `strip_version = TRUE` removes a trailing
`.<digits>` from transcript IDs (summing counts of any IDs that
collapse) to defend against this.

## See also

[`ReadSTARsoloSJ`](https://avisrilab.github.io/Matisse/reference/ReadSTARsoloSJ.md)
for the short-read (junction) path;
[`CreateMatisseObject`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md).

## Examples

``` r
if (FALSE) { # \dontrun{
txc <- ReadTranscriptMatrix("bagpiper/mats", cells = colnames(seu))
obj <- CreateMatisseObject(seu, transcript_counts = txc,
                           events = ioe_files)
} # }
```
