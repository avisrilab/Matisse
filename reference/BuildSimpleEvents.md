# Build a minimal junction event annotation table

Helper for quickly creating the `event_data` data.frame required by
[`CreateMatisseObject`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md)
and
[`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
from a table of per-junction metadata.

## Usage

``` r
BuildSimpleEvents(junctions, event_type = "simple")
```

## Arguments

- junctions:

  A `data.frame` with at minimum columns `junction_id` and `gene_id`.

- event_type:

  Character. Currently only `"simple"` is supported.

## Value

A `data.frame` suitable for the `event_data` argument of
[`CreateMatisseObject`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md).

## Details

Each row in `junctions` defines one junction. Events are then built by
pairing junctions from the same gene according to `event_type`:

- simple:

  Each junction is treated as its own event. Inclusion = the junction
  itself; exclusion = all other junctions for the same gene.

This is a convenience function for exploratory use. For production
analyses, supply a fully annotated `event_data` table (e.g. from rMATS
or SUPPA2).
