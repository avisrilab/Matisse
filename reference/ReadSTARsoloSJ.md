# Read a STARsolo splice-junction (SJ) matrix

Loads the per-cell splice-junction count matrix produced by STARsolo
with `--soloFeatures SJ` and relabels its junctions into the
`chr-start-end-strand` ID form Matisse expects for junction mode (the
first pattern recognised by the internal junction-ID parser, so the
counts are ready for
[`CreateMatisseObject`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md)
and sashimi plots).

## Usage

``` r
ReadSTARsoloSJ(
  sj_dir,
  cells = NULL,
  strand_map = c(`0` = "*", `1` = "+", `2` = "-"),
  infer_strand = TRUE,
  verbose = TRUE
)
```

## Arguments

- sj_dir:

  Directory holding `matrix.mtx`, `barcodes.tsv` and `features.tsv`
  (optionally `.gz`). Typically a STARsolo `Solo.out/SJ/raw` directory.

- cells:

  Optional character vector of cell barcodes to keep (e.g.
  `colnames(seurat)`). Barcodes absent from the SJ matrix are dropped
  with a warning. Default: `NULL` (keep all).

- strand_map:

  Named character vector mapping STAR strand codes to symbols. Default:
  `c("0" = "*", "1" = "+", "2" = "-")`.

- infer_strand:

  Logical. When the STAR strand code is `0` (undefined) but the
  intron-motif column is canonical, recover the strand from the motif
  (`1/3/5 -> +`, `2/4/6 -> -`; `0` stays undefined). Roughly one in six
  STARsolo junctions are strand-undefined; without this they can never
  match the `+`/`-` junction IDs in a SUPPA2 event table. Default:
  `TRUE`.

- verbose:

  Logical. Print progress. Default: `TRUE`.

## Value

A `dgCMatrix` of raw junction counts, *cells x junctions* (Matisse
convention), with junction IDs as column names. Feed directly to
`CreateMatisseObject(junction_counts = ...)`.

## Details

STARsolo writes the SJ matrix as *junctions x cells* with a 9-column
`features.tsv` (`chr`, intron start, intron end, strand code, splice
motif, ...). The intron start/end are the first and last intronic base;
the strand code is `0` (undefined), `1` (`+`) or `2` (`-`). Only a `raw`
matrix is emitted for SJ, so pass the filtered cell barcodes (e.g. from
`Gene/filtered`) via `cells` to subset.

## See also

[`BuildJunctionEvents`](https://avisrilab.github.io/Matisse/reference/BuildJunctionEvents.md)
to turn SUPPA2 events into a junction-ID event table that matches these
column names.

## Examples

``` r
if (FALSE) { # \dontrun{
jxn <- ReadSTARsoloSJ("Solo.out/SJ/raw", cells = colnames(seu))
obj <- CreateMatisseObject(seu, junction_counts = jxn,
                           events = BuildJunctionEvents(ioe_files))
} # }
```
